#!/usr/bin/env python3
"""
KEGG KO Annotator & Pathway Map Downloader
------------------------------------------
支持根据输入文件（TSV/CSV/Excel）中的 KO 列表检索基因名称、功能描述和 EC 编号。
支持根据指定的 KEGG Pathway Map 自动抓取对应的 KO 列表并进行注释检索。
支持直接在命令行中输入 KO 列表进行检索。
支持智能本地缓存（KO 缓存、Pathway-KO 映射缓存）。
"""

import argparse
import sys
import os
import json
import time
import re
import requests
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_KO_CACHE_FILE = os.path.join(SCRIPT_DIR, "cache", "kegg_ko_cache.json")
DEFAULT_PATHWAY_CACHE_FILE = os.path.join(SCRIPT_DIR, "cache", "kegg_pathway_cache.json")
CACHE_EXPIRE_DAYS = 7

# =====================================================================
# 1. 网络请求重试与 API 客户端 (预留扩展接口)
# =====================================================================

class KeggAPIClient:
    """KEGG REST API 客户端，带指数退避的重试机制"""
    BASE_URL = "http://rest.kegg.jp"

    def __init__(self, timeout=30, max_retries=3):
        self.timeout = timeout
        self.max_retries = max_retries

    def _get(self, endpoint):
        url = f"{self.BASE_URL}/{endpoint.lstrip('/')}"
        for attempt in range(1, self.max_retries + 1):
            try:
                resp = requests.get(url, timeout=self.timeout)
                resp.raise_for_status()
                return resp.text
            except requests.RequestException as e:
                print(f"[警告] 访问 KEGG API 失败 (尝试 {attempt}/{self.max_retries}): {e}", file=sys.stderr)
                if attempt == self.max_retries:
                    raise e
                time.sleep(2 ** attempt)  # 指数退避

    def list_ko(self):
        """获取所有 KO 定义 (KXXXXX -> definition)"""
        return self._get("list/ko")

    def list_pathway(self):
        """获取所有 Pathway 列表 (path:mapXXXXX -> name)"""
        return self._get("list/pathway")

    def link_ko_by_pathway(self, pathway_num):
        """获取指定通路下的 KO 列表"""
        target_pathway = f"path:ko{pathway_num}"
        print(f"正在从 KEGG 获取通路 {target_pathway} 的 KO 映射关系...")
        return self._get(f"link/ko/{target_pathway}")

    # 未来可在此扩展新的 API 接口，例如：
    # def list_compound(self): ...
    # def get_pathway_module(self, module_id): ...


# =====================================================================
# 2. 智能缓存管理器
# =====================================================================

class KeggCacheManager:
    """缓存管理器，用于管理 KO 数据以及 Pathway 元数据"""
    def __init__(self, ko_cache_file=DEFAULT_KO_CACHE_FILE, pathway_cache_file=DEFAULT_PATHWAY_CACHE_FILE, expire_days=CACHE_EXPIRE_DAYS):
        self.ko_cache_file = ko_cache_file
        self.pathway_cache_file = pathway_cache_file
        self.expire_days = expire_days

    def _is_cache_valid(self, cache_data):
        if not cache_data or 'timestamp' not in cache_data:
            return False
        if self.expire_days is None:
            return True
        cache_time = cache_data.get('timestamp', 0)
        return (time.time() - cache_time) < (self.expire_days * 86400)

    def load_ko_definitions(self, client, force_update=False):
        """载入 KO 定义，兼容旧版的 json 格式"""
        if not force_update and os.path.exists(self.ko_cache_file):
            try:
                with open(self.ko_cache_file, 'r', encoding='utf-8') as f:
                    cache = json.load(f)
                if self._is_cache_valid(cache):
                    print(f"[缓存] 使用本地 KO 缓存: {self.ko_cache_file}")
                    return cache['data']
                else:
                    print(f"[缓存] KO 缓存已过期，准备重新下载")
            except Exception as e:
                print(f"[警告] 读取 KO 缓存失败: {e}，将重新下载", file=sys.stderr)

        # 下载并解析
        print("正在从 KEGG API 下载所有 KO 定义（约 2 万条，几秒钟）...")
        try:
            raw_text = client.list_ko()
        except Exception as e:
            print(f"[错误] 无法下载 KO 定义: {e}", file=sys.stderr)
            sys.exit(1)

        ko_dict = {}
        for line in raw_text.strip().splitlines():
            if '\t' not in line:
                continue
            ko_id, definition = line.split('\t', 1)
            ko_id = ko_id.replace("ko:", "").strip().upper()
            ko_dict[ko_id] = definition.strip()

        # 保存到缓存
        cache_data = {
            'timestamp': time.time(),
            'data': ko_dict
        }
        try:
            with open(self.ko_cache_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, ensure_ascii=False)
            print(f"[缓存] 新的 KO 定义已保存至 {self.ko_cache_file}")
        except Exception as e:
            print(f"[警告] 保存 KO 缓存失败: {e}", file=sys.stderr)

        return ko_dict

    def load_pathway_cache(self):
        """加载 Pathway 及映射关系缓存"""
        if os.path.exists(self.pathway_cache_file):
            try:
                with open(self.pathway_cache_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception as e:
                print(f"[警告] 读取 Pathway 缓存失败: {e}", file=sys.stderr)
        return {
            "pathways": {"timestamp": 0, "data": {}},
            "pathway_ko_links": {"timestamp": 0, "data": {}}
        }

    def save_pathway_cache(self, cache_data):
        try:
            with open(self.pathway_cache_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, ensure_ascii=False, indent=2)
        except Exception as e:
            print(f"[警告] 保存 Pathway 缓存失败: {e}", file=sys.stderr)

    def get_pathways(self, client, force_update=False):
        """获取所有 Pathway (path:mapXXXXX -> name)"""
        cache = self.load_pathway_cache()
        pw_cache = cache.get("pathways", {"timestamp": 0, "data": {}})
        
        if not force_update and self._is_cache_valid(pw_cache):
            return pw_cache["data"]

        print("正在从 KEGG API 下载所有 Pathway 通路列表...")
        try:
            raw_text = client.list_pathway()
        except Exception as e:
            print(f"[错误] 无法获取 Pathway 列表: {e}", file=sys.stderr)
            return pw_cache.get("data", {})  # 失败时返回旧缓存

        pathway_dict = {}
        for line in raw_text.strip().splitlines():
            if '\t' not in line:
                continue
            pw_id, pw_name = line.split('\t', 1)
            pathway_dict[pw_id.strip()] = pw_name.strip()

        cache["pathways"] = {
            "timestamp": time.time(),
            "data": pathway_dict
        }
        self.save_pathway_cache(cache)
        return pathway_dict

    def get_pathway_kos(self, client, num_id, force_update=False):
        """获取指定通路的所有 KO 号列表"""
        target_pathway = f"path:ko{num_id}"

        cache = self.load_pathway_cache()
        links_cache = cache.setdefault("pathway_ko_links", {"timestamp": 0, "data": {}})
        pw_links = links_cache.get("data", {})
        
        # 检查该通路在缓存中是否存在且未过期
        if not force_update and target_pathway in pw_links:
            entry = pw_links[target_pathway]
            # 兼容带时间戳的子项
            if isinstance(entry, dict) and self._is_cache_valid(entry):
                print(f"[缓存] 使用本地缓存的 {target_pathway} -> KO 映射")
                return entry["data"]
            elif isinstance(entry, list):
                print(f"[缓存] 使用本地缓存的 {target_pathway} -> KO 映射")
                return entry

        # 下载
        try:
            raw_text = client.link_ko_by_pathway(num_id)
        except Exception as e:
            print(f"[错误] 无法获取 Pathway 到 KO 的映射: {e}", file=sys.stderr)
            sys.exit(1)

        ko_list = []
        for line in raw_text.strip().splitlines():
            if '\t' not in line:
                continue
            _, ko_val = line.split('\t', 1)
            ko_id = ko_val.replace("ko:", "").strip().upper()
            ko_list.append(ko_id)

        # 缓存单条映射关系，并附带时间戳
        pw_links[target_pathway] = {
            "timestamp": time.time(),
            "data": ko_list
        }
        cache["pathway_ko_links"]["data"] = pw_links
        self.save_pathway_cache(cache)
        return ko_list

# =====================================================================
# 3. 数据处理与辅助函数
# =====================================================================

def parse_ec_from_definition(definition):
    """从定义中提取 EC 号"""
    match = re.search(r'\[EC:([^\]]+)\]', definition)
    return match.group(1) if match else ""

def parse_gene_name(definition):
    """从定义中提取基因名称"""
    # 去掉 EC 标记，保留 NAME 字段完整内容
    return re.sub(r'\s*\[EC:[^\]]+\]$', '', definition).strip()

def sanitize_filename(name):
    """净化文件名，去除非法字符"""
    return re.sub(r'[\\/*?:"<>| ]', '_', name)

def read_input_file(filepath):
    """支持读取 TSV, CSV, Excel"""
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.csv':
        return pd.read_csv(filepath, dtype=str)
    elif ext in ['.xls', '.xlsx']:
        return pd.read_excel(filepath, dtype=str)
    else:
        # 默认按 TSV 读取
        return pd.read_csv(filepath, sep='\t', dtype=str)

def write_output_file(df, filepath):
    """支持写入 TSV, CSV, Excel"""
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.csv':
        df.to_csv(filepath, index=False)
    elif ext in ['.xls', '.xlsx']:
        df.to_excel(filepath, index=False)
    else:
        # 默认按 TSV 写入
        df.to_csv(filepath, sep='\t', index=False)
    print(f"数据保存完成！结果已写入 {filepath}")

# =====================================================================
# 4. 主程序入口
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="KEGG 信息抓取与 KO 注释工具 (支持 Pathway map 直接提取)")
    
    # 模式 A: 基于文件的 KO 列表注释 (Legacy)
    parser.add_argument("--input", help="输入文件路径 (支持 .tsv, .csv, .xlsx)")
    parser.add_argument("--ko_column", help="包含 KO 号的列名 (文件模式下必填)")
    
    # 模式 B: 基于 Pathway Map 的 KO 提取与注释
    parser.add_argument("--map", help="KEGG Pathway Map ID，如 map00190、ko00190 或 00190。")
    
    # 模式 C: 命令行直接检索 KO 列表
    parser.add_argument("--ko", help="直接输入的 KO 编号，逗号分隔，例如 K00001,K00002")

    # 通用选项
    parser.add_argument("--output", help="输出文件路径 (支持 .tsv, .csv, .xlsx)。如果不指定，在 Pathway 模式下将自动生成默认文件名")
    parser.add_argument("--update-cache", action="store_true", help="强制重新下载 KEGG 映射数据并更新本地缓存")
    parser.add_argument("--cache-file", default=DEFAULT_KO_CACHE_FILE, help="KO 定义缓存文件路径")
    parser.add_argument("--pathway-cache", default=DEFAULT_PATHWAY_CACHE_FILE, help="Pathway 缓存文件路径")
    parser.add_argument("--expire-days", type=int, default=CACHE_EXPIRE_DAYS, help="缓存有效期天数，设为 -1 则永不过期")
    
    args = parser.parse_args()

    # 自动创建 cache 目录
    os.makedirs(os.path.dirname(args.cache_file), exist_ok=True)
    os.makedirs(os.path.dirname(args.pathway_cache), exist_ok=True)
    os.makedirs(os.path.join(SCRIPT_DIR, "data", "output"), exist_ok=True)
    
    # 规范化缓存过期设置
    expire_days = None if args.expire_days == -1 else args.expire_days

    # 初始化客户端与缓存管理器
    client = KeggAPIClient()
    cache_mgr = KeggCacheManager(
        ko_cache_file=args.cache_file, 
        pathway_cache_file=args.pathway_cache, 
        expire_days=expire_days
    )

    # 载入所有 KO 映射表定义 (约 2w 条数据，支持本地缓存)
    ko_def_map = cache_mgr.load_ko_definitions(client, force_update=args.update_cache)

    df = None
    output_path = args.output

    # ================= 流程判断 =================
    if args.map:
        # --- 模式 B: Pathway Map 检索 ---
        pw_id = args.map.strip()
        
        # 1. 载入 pathway 列表以获取通路中文/英文描述名称
        pathway_names = cache_mgr.get_pathways(client, force_update=args.update_cache)
        
        # 2. 解析出标准的 map 数字 ID 键名
        match = re.search(r'\d{5}', pw_id)
        if not match:
            print(f"[错误] 无效的 Pathway Map ID '{pw_id}'，必须包含 5 位数字 (例如 map00010 或 00010)。", file=sys.stderr)
            sys.exit(1)
        num_id = match.group(0)
        
        # 寻找对应的标准名
        std_pw_id = f"map{num_id}"
        pw_name = pathway_names.get(std_pw_id, pathway_names.get(f"path:{std_pw_id}", f"Pathway_{num_id}"))
        print(f"匹配到通路: {std_pw_id} - {pw_name}")

        # 3. 获取该 Pathway 对应的 KO 号列表
        ko_list = cache_mgr.get_pathway_kos(client, num_id, force_update=args.update_cache)
        print(f"通路 {std_pw_id} 共有 {len(ko_list)} 个 KO 成员。")

        if not ko_list:
            print(f"[警告] 通路 {std_pw_id} 没有关联的 KO。可能是没有该通路或 KEGG 暂无数据。")
            df = pd.DataFrame(columns=['KO'])
        else:
            df = pd.DataFrame({'KO': ko_list})

        # 4. 自动推断默认输出文件名 (如果用户没有提供 --output)
        if not output_path:
            safe_name = sanitize_filename(pw_name)
            output_path = os.path.join(SCRIPT_DIR, "data", "output", f"{safe_name}.tsv")
            print(f"未指定输出文件，自动命名并保存至: {output_path}")

    elif args.ko:
        # --- 模式 C: 命令行直接指定 KO ---
        ko_list = [k.strip().upper() for k in args.ko.split(',') if k.strip()]
        print(f"从命令行获取到 {len(ko_list)} 个待查询 KO。")
        df = pd.DataFrame({'KO': ko_list})
        
        if not output_path:
            output_path = os.path.join(SCRIPT_DIR, "data", "output", "ko_annotations.tsv")
            print(f"未指定输出文件，默认保存为: {output_path}")

    elif args.input:
        # --- 模式 A: 传统文件输入模式 ---
        if not args.ko_column:
            print("[错误] 文件模式下，必须指定包含 KO 编号的列名 (--ko_column)", file=sys.stderr)
            sys.exit(1)
        
        try:
            df = read_input_file(args.input)
            print(f"成功读取输入文件 {args.input}，共 {len(df)} 行。")
        except Exception as e:
            print(f"[错误] 无法读取输入文件 {args.input}: {e}", file=sys.stderr)
            sys.exit(1)

        if args.ko_column not in df.columns:
            print(f"[错误] 列名 '{args.ko_column}' 不在输入文件中。可用列有: {list(df.columns)}", file=sys.stderr)
            sys.exit(1)

        if not output_path:
            base, ext = os.path.splitext(args.input)
            # 如果输入在 data/input/，我们输出到 data/output/
            input_dir, input_name = os.path.split(args.input)
            if os.path.basename(input_dir) == "input" and os.path.basename(os.path.dirname(input_dir)) == "data":
                output_path = os.path.join(SCRIPT_DIR, "data", "output", f"{os.path.splitext(input_name)[0]}_annotated{ext}")
            else:
                output_path = f"{base}_annotated{ext}"
            print(f"未指定输出文件，默认命名为: {output_path}")
    else:
        parser.print_help()
        sys.exit(0)

    # ================= 统一检索注释逻辑 =================
    if df is not None:
        ko_col = args.ko_column if (args.input and args.ko_column) else 'KO'
        
        definitions = []
        gene_names = []
        ec_numbers = []

        for ko in df[ko_col]:
            if pd.isna(ko):
                definitions.append("")
                gene_names.append("")
                ec_numbers.append("")
                continue

            ko_clean = str(ko).strip().upper()
            if ko_clean.startswith("KO:"):
                ko_clean = ko_clean[3:]
            
            definition = ko_def_map.get(ko_clean, "")
            definitions.append(definition)
            if definition:
                gene_names.append(parse_gene_name(definition))
                ec_numbers.append(parse_ec_from_definition(definition))
            else:
                gene_names.append("")
                ec_numbers.append("")

        df['definition'] = definitions
        df['gene_name'] = gene_names
        df['ec_number'] = ec_numbers

        # 保存结果
        write_output_file(df, output_path)

if __name__ == "__main__":
    main()
