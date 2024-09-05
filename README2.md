
#プログラムをダウンロードしてからすること。
1, ダウンロードしたものがあるディレクトリに移動します。
2, bash ch_shebang.sh を実行します。
3, chmod +x *.py を実行します。
4, ku1mV.py init を実行します。

#解析手順 ( 仮に240905 に観測したSN_2024cld を解析することを想定 )
1, cd ${HOME}/ku1mV/kSIRIUS を実行して
    mkdir 240905 と cd 240905 を実行して、
    ここに 240905 の kSIRIUS のデータをダウンロード(scp のやつ)
    cd ${HOME}/ku1mV/optcam を実行して
    mkdir 240905 と cd 240905 を実行して、
    ここに 240905 の 可視カメラ のデータをダウンロード (scp のやつ)
2, ダウンロードしたものがあるディレクトリに移動します。
3, ku1mV.py <天体名> を実行します。ex) ku1mV.py SN_2024cld
4, ku1mV.py <YYMMDD> <天体名> で処理を開始します。ex) ku1mV.py 240905 SN_2024cld