
0, 何もわからない人はこっち [README2.md の内容を見る](README2.md)

1, まず、ch_shebang.sh を実行します。
2, main.param を適切に設定します。
3, main.py <任意のobject名> を実行します。
4, '任意のobject名' ファイルを適切に編集します。
5, advanced.param をいい感じに編集します。ほとんどの場合不要。
6, ku1mV.py <YYMMDD>  <'任意のobject名'> で処理を開始します。

詳細：
1
    ch_shebang.sh を行うことによって、すべてのプログラムのシバンを適切に変更します。
    これは、python3 がどこにインストールされているかを示すものです。
    仮想環境等は active にした状態で実行してください。

2  
    main.param の中では、5 で使用する、
    YYMMDD, '任意のobject名' がそれぞれ、{date}, {objectname} として使用可能です。
    objfile_dir は '任意のobject名' ファイルを作成するディレクトリを指定します。
    rawdata_dir は 未処理の画像が保存されているディレクトリを指定します。
    work_dir は 処理が行われるディレクトリを指定します。
    quicklook は実装されていません。
    flatdiv はフラットで割る作業をします。

3
    '任意のobject名' はどんな名前でも構いません。
    天体ごとにまとめて作業するための名前づけです。

4
    '任意のobject名' ファイル は main.param で指定した。
    objfile_dir に保存されているはずです。ここで、fits 画像の header 
    に登録されている OBJECT の、一部(他天体と区別可能な範囲)を追記することで、
    適切に fits 画像を検索できます。
    いくつか追加することで観測者のミス等による表記揺れの対策ができます。

5
   一番気になっているのは maxstarnum と minstarnum です。
   それぞれ、match の計算に時間がかかる時 と match できない時 に変えます。
   
6 
   YYYYMMDD でも行けるといいな