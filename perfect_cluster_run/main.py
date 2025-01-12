import subprocess
import re
import os
from concurrent.futures import ProcessPoolExecutor
import sys

repeat = sys.argv[1]
problem = sys.argv[2]
# ell = [36, 72, 108, 144]
# ell = [204, 216, 228]
# ell = [84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204]

ell = [36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228]

print(f"problem id = {problem}")

# 0 onemax
# 1 mk
# 2 FTRAP
# 3 cyc
# 4 nk
# 5 SPIN
# 6 SAT
# 7 3-opt

# 設定參數範圍


# ell = [9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153, 156, 159, 162, 165, 168, 171, 174, 177, 180]  # 可擴展到其他值，如 48, 60, ...
# ell = [9, 12, 15, 18, 21, 24, 27, 30]  # 可擴展到其他值，如 48, 60, ...
# ell = [33, 36, 39, 42, 45, 48, 51, 54, 57, 60]  # 可擴展到其他值，如 48, 60, ...
# ell = [63, 66, 69, 72, 75, 78, 81, 84, 87, 90]  # 可擴展到其他值，如 48, 60, ...
# ell = [93, 96, 99, 102, 105, 108, 111, 114, 117, 120]  # 可擴展到其他值，如 48, 60, ...
# ell = [123, 126, 129, 132, 135, 138, 141, 144, 147, 150]  # 可擴展到其他值，如 48, 60, ...
# ell = [153, 156, 159, 162, 165, 168, 171, 174, 177, 180]  # 可擴展到其他值，如 48, 60, ...
# ell = [12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204]  # 可擴展到其他值，如 48, 60, ...

# repeat = repeat ## check this!
# problem = problem



# 結果存放資料夾
output_folder = "0_folder"
os.makedirs(output_folder, exist_ok=True)

# 存放結果的字典
results = {}

def run_command(ell):
    """執行指令，提取 NFE，並將結果存檔"""
    if problem == 6:
        command = f"./sweep {ell} {repeat} {problem} 0"
    else:
        command = f"./sweep {ell} {repeat} {problem}"

    try:
        # 執行指令並捕獲輸出
        output = subprocess.check_output(command, shell=True, text=True)

        # 使用正則表達式提取 NFE
        match = re.search(r"NFE:\s*([\d.]+)", output)
        if match:
            nfe = float(match.group(1))
            results[ell] = nfe

            # 儲存結果到對應的檔案
            output_file = os.path.join(output_folder, f"{ell}.txt")
            with open(output_file, "w") as f:
                f.write(f"ell: {ell}, NFE: {nfe}\n")

            print(f"ell: {ell}, NFE: {nfe}")
        else:
            print(f"未找到 NFE 值: {output}")
    except subprocess.CalledProcessError as e:
        print(f"執行失敗: {e}")

# 使用 ProcessPoolExecutor 進行平行執行
with ProcessPoolExecutor() as executor:
    executor.map(run_command, ell)
