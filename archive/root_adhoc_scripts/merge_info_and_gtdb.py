import csv
import os

# 定义输入输出路径（使用绝对路径以防在不同目录运行脚本时找不到文件）
project_root = "/dssg/home/acct-schhui/gses0/rentianxin"
info_file = os.path.join(project_root, "data/gses_all_rawMAGs_fa_V0/all_success_sample_MAG_info.tsv")
gtdb_file = os.path.join(project_root, "results/phylogeny/gtdbtk/gses_all_mags_V0/merged_gtdbtk_summary.tsv")
output_file = os.path.join(project_root, "results/phylogeny/gtdbtk/gses_all_mags_V0/merged_mag_info_gtdbtk.tsv")

print("正在读取 GTDB-Tk 分类汇总文件...")
gtdb_data = {}
with open(gtdb_file, "r", encoding="utf-8") as f:
    # 使用 DictReader 自动识别列名并以 tab 分隔
    reader = csv.DictReader(f, delimiter="\t")
    gtdb_fields = reader.fieldnames
    for row in reader:
        # 以 user_genome 列作为唯一主键
        gtdb_data[row["user_genome"]] = row

print("正在对齐并拼接 MAG 信息...")
merged_rows = []
merged_fields = []

with open(info_file, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    info_fields = reader.fieldnames
    
    # 构造合并后的表头（保留 MAG info 的列，后面拼接 GTDB-Tk 除了 user_genome 之外的所有列）
    gtdb_fields_subset = [field for field in gtdb_fields if field != "user_genome"]
    merged_fields = info_fields + gtdb_fields_subset
    
    for row in reader:
        bin_id = row["bin"]
        gtdb_row = gtdb_data.get(bin_id)
        
        merged_row = {}
        # 先复制已成功样本的 MAG info
        merged_row.update(row)
        
        if gtdb_row:
            # 找到对应的 GTDB-Tk 分类信息，横向拼入
            for field in gtdb_fields_subset:
                merged_row[field] = gtdb_row.get(field, "N/A")
        else:
            # 如果没有找到对应的分类信息，则填充为 N/A 避免数据缺失
            for field in gtdb_fields_subset:
                merged_row[field] = "N/A"
                
        merged_rows.append(merged_row)

print(f"正在写入合并后的文件: {output_file} ...")
# 确保输出目录存在
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=merged_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(merged_rows)

print(f"完成！成功合并了 {len(merged_rows)} 行数据。")
