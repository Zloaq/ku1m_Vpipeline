#!/opt/anaconda3/envs/p11/bin/python3

import subprocess
import re
import os

def save_directory_names(file_path):
    dir_names = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # 対象となるディレクトリ名のパターン
            if re.match(r'^(objfile_dir|rawdata_infra|rawdata_opt|work_dir)', line):
                key, value = line.split(maxsplit=1)
                dir_names[key] = value.strip()  # keyと値を保存
    return dir_names


def git_pull():
    try:
        print("Pulling the latest updates from the repository...")
        result = subprocess.run(['git', 'reset', '--hard'], check=True, text=True, capture_output=True)
        print(result.stdout)  # git pull の結果を表示
        result = subprocess.run(['git', 'pull'], check=True, text=True, capture_output=True)
        print(result.stdout)  # git pull の結果を表示
    except subprocess.CalledProcessError as e:
        print(f"Error during git pull: {e.stderr}")


def restore_directory_names(file_path, dir_names):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            key = line.split(maxsplit=1)[0] if line.strip() else ""
            if key in dir_names:
                file.write(f"{key:<16}{dir_names[key]}\n")
            else:
                file.write(line)


if __name__ == "__main__":

    path_program = os.path.abspath(__file__)
    dir_of_program = os.path.dirname(path_program)
    param_file = os.path.join(dir_of_program, 'main.param')

    saved_dir_names = save_directory_names(param_file)

    git_pull()

    restore_directory_names(param_file, saved_dir_names)

    print("Update completed and directory names restored.")