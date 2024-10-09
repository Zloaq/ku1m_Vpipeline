

0, 何もわからない人はこっち [README2.md の内容を見る](README2.md)

1, まず、git clone https://github.com/Zloaq/ku1m_Vpipeline を実行する。
2, ch_shebang.sh を実行します。
3, main.param を適切に設定します。
4, ku1mV.py <任意のobject名> を実行します。
5, '任意のobject名' ファイルを適切に編集します。
6, advanced.param をいい感じに編集します。ほとんどの場合不要。
7, ku1mV.py <YYMMDD>  <'任意のobject名'> で処理を開始します。

・ku1m init でホームディレクトリにディレクトリ構造が作成され、main.param もそれに合わせて変更されます。
　作成されるディレクトリ名は ku1mV で中身は
　生画像ディレクトリ、作業ディレクトリ、object名ファイルディレクトリです。

・update_pipeline.py を実行すると、github上のアップデートを適応できます。
　main.param のディレクトリ名は維持されます。


詳細：
1
    git clone https://github.com/Zloaq/ku1m_Vpipeline を実行した場所に
    ku1m_Vpipeline というディレクトリが作成され、
    その中にリポジトリの内容がダウンロードされます。
    git clone https://github.com/Zloaq/ku1m_Vpipeline <ディレクトリ名>
    で好きなディレクトリ名でダウンロードできます。
    git がない人は sudo apt install git
    

2
    ch_shebang.sh を行うことによって、すべてのプログラムのシバンを適切に変更します。
    これは、python3 がどこにインストールされているかを示すものです。
    仮想環境等は active にした状態で実行してください。

3  
    main.param の中では、5 で使用する、
    YYMMDD, '任意のobject名' がそれぞれ、{date}, {objectname} として使用可能です。
    objfile_dir は '任意のobject名' ファイルを作成するディレクトリを指定します。
    rawdata_dir は 未処理の画像が保存されているディレクトリを指定します。
    work_dir は 処理が行われるディレクトリを指定します。
    quicklook は実装されていません。
    flatdiv はフラットで割る作業をします。

4
    '任意のobject名' はどんな名前でも構いません。
    天体ごとにまとめて作業するための名前づけです。

5
    '任意のobject名' ファイル は main.param で指定した。
    objfile_dir に保存されているはずです。ここで、fits 画像の header 
    に登録されている OBJECT の、一部(他天体と区別可能な範囲)を追記することで、
    適切に fits 画像を検索できます。
    いくつか追加することで観測者のミス等による表記揺れの対策ができます。

6
   一番気になっているのは maxstarnum と minstarnum です。
   それぞれ、match の計算に時間がかかる時 と match できない時 に変えます。
   
7 
   YYYYMMDD でも行けるといいな