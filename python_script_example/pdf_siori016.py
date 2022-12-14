import os,sys
from tkinter import *
from tkinter import ttk
from tkinter import filedialog

"""
01.PDFから栞のリストを取出す 220401
02.取出したリスト通りに栞を挿入しなおす。220401
03.PDFからPDFへ01をした後02を実行。221110
以下参考サイト
"""
#Tkinterに表示したタブの中に、さらにタブを作成したい。
#https://teratail.com/questions/rxqmt75cfd3lme

#Tkinterの階層構造とフレーム（Frame）
#https://denno-sekai.com/tkinter-frame/

#NOTE t.koba(2021.3.9)
#https://note.com/uranus_xii_jp/n/nbf8c42c0625b

#以下PDFの栞の書き込み
#https://jablogs.com/detail/49869

import fitz
import csv
from PyPDF2 import PdfFileWriter, PdfFileReader

def pdf_to_csv(pdf_in,csv_out):
    doc = fitz.open(pdf_in)
    page_count = doc.page_count # 総ページ数
    print('Page count: ', page_count)
    toc_list = doc.get_toc()

    with open(csv_out, 'w',newline="") as f: #空白行をなくすためnewline追加(220329)
        writer = csv.writer(f)
        writer.writerows(toc_list) #org
    f.close()
    return()

def csv_to_pdf(pdf_in, csv_out):
    with open(csv_out) as f: #csvの読込み
        reader = csv.reader(f, delimiter=',')
        l = [row for row in reader]
    
    l2 = [] #csvの空データ行削除
    for i in range(len(l)):
        if l[i] == []:
            pass
        else:
            #print(l[i])
            l2.append(l[i])

    output = PdfFileWriter() # open output
    input1 = PdfFileReader(open(pdf_in, 'rb')) # open input
    
    n = input1.getNumPages()
    for i in range(n):
        output.addPage(input1.getPage(i)) # insert page

    par_list = [] #階層のある場合のための階層親リスト
    kaiso_list = [] #階層番号のリスト

    for i in l2:
        kaiso = int(i[0])
        kaiso_list.append(kaiso)
        name = i[1]
        page = int(i[2])
        #print(kaiso,name,page)
        #page は -1 しないとずれる(2022.03.28)
        if kaiso == 1:
            par_list.append(output.addBookmark(name, page-1, parent=None))
        else:
            ln = len(kaiso_list)
            par_kaiso = kaiso-1 #親の階層
            r_kaiso = kaiso_list[::-1] #階層リストの逆順
            ind = r_kaiso.index(par_kaiso) #逆から親の階層検索
            p_ind = ln - ind - 1 #par_listの親階層インデクス
            par = par_list[p_ind]
            par_list.append(output.addBookmark(name, page-1, parent=par))

    base_pdf_in = os.path.splitext(os.path.basename(pdf_in))[0] #拡張子なしのファイル名
    sname = base_pdf_in + "_栞.pdf"
    #outputStream = open('result.pdf','wb') #creating result pdf JCT
    outputStream = open(sname,'wb') #221105
    output.write(outputStream) #writing to result pdf JCT
    outputStream.close() #closing result JCT

#221108開発 PDF[0]の栞をPDF[1]に適用する。
def pdf_to_pdf(pdf_in0, pdf_in1):
    csv_in0 = os.path.splitext(os.path.basename(pdf_in0))[0] + ".csv" #拡張子なしのファイル名
    pdf_to_csv(pdf_in0,csv_in0)
    print(pdf_in0)
    csv_to_pdf(pdf_in1,csv_in0)
    print("完了 pdf[0]->pdf[1]")

class PdfGui():
    def __init__(self):
        dirFile = os.path.abspath(os.path.dirname(__file__)) #ファイルのある場所
        os.chdir(dirFile) #ファイルのある場所にcd移動
        os.chdir('../') #カレントディレクトリを1つ上に移動する。
        self.iFile = os.getcwd() #カレントディレクトリのパス保管
        root = Tk() # rootの作成
        self.root = root
        root.title("PDFtoCSV16 kitagawa")
        root.geometry("500x220") #メインウィンドウの大きさを設定
        nb = ttk.Notebook(root) #メインウィンドウにnotebookを作成する。(タブの作成)
        tab1 = Frame(nb) #notebookに関するフレームを作る。
        tab2 = Frame(nb)
        nb.add(tab1, text="PDF-CSV", padding=3) #notebookに対してtab1, 2をそれぞれ追加する。
        nb.add(tab2, text="PDF-PDF", padding=3)
        nb.pack(expand=1, fill="both") #メインフレームでのnotebook配置を決定する。
        self.tab1_main(tab1) #各タブの内容を記載する。
        self.tab2_main(tab2)
        root.mainloop()

    def tab1_main(self, tab1):
        nb1 = tab1
        frame1 = ttk.Frame(nb1, padding=10)
        frame1.grid(row=1, column=1, sticky=E)
        IFileLabel = ttk.Label(frame1, text="ファイル参照(pdf)＞＞", padding=(5, 2)) #「ファイル参照」ラベルの作成
        IFileLabel.pack(side=LEFT)
        entry1 = StringVar() # 「ファイル参照」エントリーの作成
        IFileEntry = ttk.Entry(frame1, textvariable=entry1, width=30)
        self.entry1 = entry1
        IFileEntry.pack(side=LEFT)
        IFileButton = ttk.Button(frame1, text="参照pdf", command=self.filedialog_clicked) #「ファイル参照」ボタンの作成
        IFileButton.pack(side=LEFT)

        # Frame2の作成
        frame2 = ttk.Frame(nb1, padding=10)
        frame2.grid(row=2, column=1, sticky=E)
        IFileLabel = ttk.Label(frame2, text="変換＞＞", padding=(5, 2))
        IFileLabel.pack(side=LEFT)
        style = ttk.Style()
        style.theme_use("winnative")
        style.configure("office.TCombobox", selectbackground="blue", padding=5)
        module = ("pdf->csv", "csv->pdf")
        v = StringVar()
        self.combobox = ttk.Combobox(frame2,textvariable= v, values=module, style="office.TCombobox")
        self.combobox.bind('<<ComboboxSelected>>', self.select_combo)
        self.combobox.pack()
        self.combobox.set("csv->pdf") 
        
        frame3 = ttk.Frame(nb1, padding=10) # Frame6の作成
        frame3.grid(row=7, column=1, sticky=E)
        IFileLabel = ttk.Label(frame3, text="ファイル参照(csv)＞＞", padding=(5, 2)) # 「ファイル参照」ラベルの作成
        IFileLabel.pack(side=LEFT)
        # 「ファイル参照」エントリーの作成
        entry3 = StringVar()
        IFileEntry = ttk.Entry(frame3, textvariable=entry3, width=30)
        self.entry3 = entry3
        IFileEntry.pack(side=LEFT)
        # 「ファイル参照」ボタンの作成
        IFileButton = ttk.Button(frame3, text="参照csv", command=self.filedialog_clicked2)
        IFileButton.pack(side=LEFT)
        
        frame4 = ttk.Frame(nb1, padding=10) # Frame7の作成
        frame4.grid(row=9,column=1,sticky=W)
        button1 = ttk.Button(frame4, text="実行", command=self.conductMain) # 実行ボタンの設置
        button1.pack(fill = "x", padx=30, side = "left")
        button2 = ttk.Button(frame4, text=("閉じる"), command=self.quit) # キャンセルボタンの設置
        button2.pack(fill = "x", padx=30, side = "left")

    def tab2_main(self, tab1): #PDFからPDFへの付箋の移動（中間ファイルとしてcsvが作成される）
        nb1 = tab1
        frame1 = ttk.Frame(nb1, padding=10)
        frame1.grid(row=1, column=1, sticky=E)
        IFileLabel = ttk.Label(frame1, text="ファイル参照(pdf[0])＞＞", padding=(5, 2)) # 「ファイル参照」ラベルの作成
        IFileLabel.pack(side=LEFT)
        entry1 = StringVar() # 「ファイル参照」エントリーの作成
        IFileEntry = ttk.Entry(frame1, textvariable=entry1, width=30)
        self.entry11 = entry1
        IFileEntry.pack(side=LEFT)
        IFileButton = ttk.Button(frame1, text="参照pdf", command=self.filedialog_clicked3) # 「ファイル参照」ボタンの作成
        IFileButton.pack(side=LEFT)

        frame2 = ttk.Frame(nb1, padding=10) # Frame6の作成
        frame2.grid(row=7, column=1, sticky=E)
        IFileLabel = ttk.Label(frame2, text="ファイル参照(pdf[1])＞＞", padding=(5, 2)) # 「ファイル参照」ラベルの作成
        IFileLabel.pack(side=LEFT)
        entry3 = StringVar() # 「ファイル参照」エントリーの作成
        IFileEntry = ttk.Entry(frame2, textvariable=entry3, width=30)
        self.entry31 = entry3
        IFileEntry.pack(side=LEFT)
        IFileButton = ttk.Button(frame2, text="参照pdf", command=self.filedialog_clicked4) # 「ファイル参照」ボタンの作成
        IFileButton.pack(side=LEFT)

        frame3 = ttk.Frame(nb1, padding=10) # Frame7の作成
        frame3.grid(row=9,column=1,sticky=W)
        button1 = ttk.Button(frame3, text="実行", command=self.conductMain2) # 実行ボタンの設置
        button1.pack(fill = "x", padx=30, side = "left")
        button2 = ttk.Button(frame3, text=("閉じる"), command=self.quit) # キャンセルボタンの設置
        button2.pack(fill = "x", padx=30, side = "left")

    def filedialog_clicked(self): #INPUT FILE
        fTyp = [("", "*")]
        iFilePath = filedialog.askopenfilename(filetype = fTyp, initialdir = self.iFile)
        self.entry1.set(iFilePath)

    def filedialog_clicked2(self): #POST FILE
        fTyp = [("", "*")]
        iFilePath = filedialog.askopenfilename(filetype = fTyp, initialdir = self.iFile)
        self.entry3.set(iFilePath)

    def filedialog_clicked3(self): #INPUT FILE
        fTyp = [("", "*")]
        iFilePath = filedialog.askopenfilename(filetype = fTyp, initialdir = self.iFile)
        self.entry11.set(iFilePath)

    def filedialog_clicked4(self): #POST FILE
        fTyp = [("", "*")]
        iFilePath = filedialog.askopenfilename(filetype = fTyp, initialdir = self.iFile)
        self.entry31.set(iFilePath)

    def conductMain(self): # 実行ボタン押下時の実行関数（タブ１）
        filePath1 = self.entry1.get() #INPUT FILE
        filePath2 = self.entry3.get() #POST FILE
        method = self.combobox.get() #変換方法
        print(filePath1)
        print(filePath2)
        self.GUI_main(filePath1,filePath2,method)

    def conductMain2(self): # 実行ボタン押下時の実行関数（タブ２）
        filePath1 = self.entry11.get() #INPUT FILE
        filePath2 = self.entry31.get() #POST FILE
        print(filePath1)
        print(filePath2)
        #self.GUI_main2(filePath1,filePath2)
        pdf_to_pdf(filePath1, filePath2)
        print("実行完了","PDF->PDF")

    def select_combo(self, event):
        print(self.combobox.get())

    def quit(self): #ウィンドウの完全な削除
        self.root.quit()
        sys.exit()

    def GUI_main(self, file_pdf,file_csv,method): #タブ１の関数
        base_file_pdf = os.path.splitext(os.path.basename(file_pdf))[0] #拡張子なしのファイル名
        file_csv0 = base_file_pdf + ".csv" #拡張子なしのファイル名
        if file_csv == "":
            file_csv = file_csv0
        else:
            None
            
        if method == "pdf->csv":
            pdf_to_csv(file_pdf,file_csv) #01. PDFから付箋データをcsvに変換する場合に使用
        elif method == "csv->pdf":
            csv_to_pdf(file_pdf,file_csv) #02. 221105 csv付箋データをPDFに書き込む時使用。書き込んだファイルは元ファイル名+"_栞".pdfで保存とした。
        else:
            print("なにもしない")
        print("実行完了",method)

PdfGui()

#if __name__ == "__main__":
#    main()