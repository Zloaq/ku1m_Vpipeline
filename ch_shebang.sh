#!/bin/bash

# 使用法: ./update_shebang.sh path_to_python_script



PYTHON_PATH=$(which python3)

SCRIPT_PATH=($(ls |grep py\$))

for filename in ${SCRIPT_PATH[@]}; do

    # シバン行を確認して更新する
    if [ -f "$filename" ]; then
        # シバン行を取得し、Pythonパスが含まれているかチェック
        SHEBANG=$(head -n 1 "$filename")
        if [[ "$SHEBANG" =~ ^#!.*python ]]; then
            # ファイルの残りの内容を保持
            tail -n +2 "$filename" > /tmp/tmpfile
            # 新しいシバン行をファイルに挿入
            echo "#!$PYTHON_PATH" > "$filename"
            cat /tmp/tmpfile >> "$filename"
            rm /tmp/tmpfile
            echo "Shebang updated to '$PYTHON_PATH' in '$filename'"
        else
            echo "No Python shebang line to update in '$filename'."
        fi
    else
        echo "File does not exist: $filename"
        exit 1
    fi

done