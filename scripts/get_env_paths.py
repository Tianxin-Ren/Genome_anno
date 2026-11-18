import subprocess
import os
import re

def main():
    print("#" + "="*70)
    print("# Snakemake 9.x 环境路径自动捕获工具 (核心版)")
    print("#" + "="*70 + "\n")

    # 运行 snakemake 的 dry-run 模式来获取环境信息
    # 我们使用 --conda-create-envs-only 确保它只关注环境
    cmd = ["snakemake", "--workflow-profile", "workflow/profiles/slurm", "--conda-create-envs-only", "--dry-run"]
    
    print(f"正在执行: {' '.join(cmd)} ...")
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        output = result.stdout
        
        # 在输出中寻找类似 .conda/163bdfb..._ 的字符串
        found_envs = re.findall(r"(\.conda/[0-9a-f]{32}_)", output)
        # 去重
        found_envs = sorted(list(set(found_envs)))

        if not found_envs:
            print("未能从输出中提取到环境路径。请确保 workflow/envs/ 目录下的 YAML 文件被 smk 引用了。")
            print("Snakemake 输出如下：")
            print(output)
            return

        print("\n" + "#"*30)
        print("# 检测到以下环境安装路径：")
        print("#"*30 + "\n")

        # 尝试匹配 YAML 文件名
        # 因为 dry-run 输出中通常会有 'Conda environment workflow/envs/xxx.yaml'
        for env_path in found_envs:
            # 这是一个简单的查找匹配，如果输出中有对应的名称则打印出来
            pattern = rf"Conda environment (workflow/envs/[^:]+):.*{env_path}"
            match = re.search(pattern, output, re.DOTALL)
            
            if match:
                yaml_file = match.group(1)
                print(f"### 环境: {os.path.basename(yaml_file)}")
                print(f"conda env create --file {yaml_file} --prefix {env_path}")
            else:
                print(f"### 发现哈希环境 (不确定对应哪个 YAML):")
                print(f"conda env create --file <对应的YAML> --prefix {env_path}")
            print("-" * 40)

    except Exception as e:
        print(f"执行出错: {e}")

if __name__ == "__main__":
    main()
