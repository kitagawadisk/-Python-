# -*- coding: utf-8 -*-
from shutil import ReadError
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import japanize_matplotlib
import os #保存先のディレクトリ作成
import csv
import re
import datetime
now = datetime.datetime.now()
now = now.strftime("%y%m%d_%H%M")
plot_color_list = ["b", "g", "r", "c", "m", "y", "k", "w"]

"""
verup情報 ver7
    不要スクリプトの削除
    importの最小化
    入力の単純化
    保存先のフォルダ作成

ver8 220518
    地山のコンタ確認

ver9
    地山レベルとボーリング支持地盤レベルの比較を可能にしたい。
    グラフの描画と計算の分離

ver10
    コンタ図のテキストより線を勝ちとすること可能に
    make_z内にある各位置補間データ作成を独立化
    post_data.csvについてgrid列を追加しon,offの読み取りを行う。それにより建物ごとのTrue,Falseリストを削除
    予定
    　secデータの作成のため、sec.csvから読込ができるようにする。
    　〇もしくは地山のデータについてもテキスト表示できるようにしておく。
    　〇各位置補間データの値については、どこを基準にするか変更できるようにする（現状20.6としているが）
    　　csvに基準を入力するようにした。
    　〇断面作成用のcsv読み込み関数の作成
    　点の分割ピッチを算定して求めるか、ピッチから分割をするか考える。
    　postファイル.csvのリスト化したものを置く。
"""

def make_sokusen(list2, n): #2次元平面の直線を線形補間で作成するだけの関数
    #測線のxy座標データを用意するためだけに作成した。220518
    #入力list2 = [[x0, y0], [x1, y1], bool]
    list22 = list2[:2] #bool排除[[p0], [p1], bool]
    list22 = list(zip(*list22)) #転置 [[x0, x1], [y0, y1]]
    xl = list22[0]
    yl = list22[1]
    y_latent = np.linspace(min(yl), max(yl), n)
    fitted_curve = interpolate.interp1d(yl, xl)
    x_latent = fitted_curve(y_latent)
    #print(x_latent)
    return(x_latent, y_latent)

def ArrayAppend(x,xi): #一次元+一次元
    x = list(x)
    for i in range(len(xi)):
        x.append(xi[i])
    x = np.array(x)
    return x

class sokusen(): #測線の計算のためだけの関数
    def __init__(self, data2, key):
        s = self
        s.key = key
        s.x0 = data2[0][0]
        s.y0 = data2[0][1]
        s.x1 = data2[1][0]
        s.y1 = data2[1][1]
        if s.key:
            s.a = ((s.y1 - s.y0)/(s.x1 - s.x0))
        else:
            s.a = ((s.x1 - s.x0)/(s.y1 - s.y0))
    
    def get_y(self, x):
        s = self
        if s.key:
            return((s.y1 - s.y0)/(s.x1 - s.x0)*(x - s.x0) + s.y0)
        else:    
            return((s.x1 - s.x0)/(s.y1 - s.y0)*(x - s.y0) + s.x0)

def kill_data(csv_file, s): #3列目のデータはレイヤ名。特定レイヤ名以外のデータは排除
    df = pd.read_csv(csv_file, header=None)
    x1 = df[0].values # = L
    z1 = df[1].values # 対応するz高さ
    k1 = df[2].values # レイヤー名
    x2 = []
    z2 = []
    k2 = []
    for i in range(len(k1)):
        #print(s, k1[i])
        if str(k1[i]) == s:
            x2.append(x1[i])
            z2.append(z1[i])
            k2.append(k1[i])
    return(x2, z2, k2)

def sokusen_xyzdata(lname, RSD):
    #測線のdxfファイルをデータ化したcsvファイルがある。
    #このCSVファイルは x方向, z方向, レイヤ名　の列名となっている。
    #支持地盤のレイヤ名は13としているため、13以外のデータを排除する。
    #また、x方向は実際は始点から終点を結ぶ測線の向かう方向の線であるため、
    #座標変換を行っている。
    #print(s_dictlist)
    #測線の支点を基準にして測線の傾き方向にL進んだ場合の
    #x,y座標を取得する
    name = lname
    s_dicti = RSD.s_dict[lname]
    s_dkeyi = RSD.s_dict[lname][2]
    #print(s_dicti)
    sk = sokusen(s_dicti, s_dkeyi) #Trueは通常関数 Falseは90度に近い線分
    #繰返し用部分
    csv_file = './csv/' + name + '.csv'
    d = kill_data(csv_file, "13") #データの中でもどのレイヤの値を使用するか。コメント時点では13レイヤのデータのみ使用。
    x2 = []
    y2 = []
    z2 = d[1]
    for Li in d[0]:
        if s_dkeyi:
            x2i = Li*np.sqrt(1./(1.+ sk.a**2)) + sk.x0
            y2i = sk.get_y(x2i)
            x2.append(x2i)
            y2.append(y2i)
        else:
            y2i = Li*np.sqrt(1./(1.+ sk.a**2)) + sk.y0
            x2i = sk.get_y(y2i)
            x2.append(x2i)
            y2.append(y2i)
    #print(len(x2),len(y2),len(z2))
    return(x2, y2, z2, d[0])

def sokusen_xyz_all():
    x2 = []
    y2 = []
    z2 = []
    RSD = Read_Section_Data("./csv/sec_測線.csv") #測線断面線データより測線ファイルのタイトルのリストを取出したい。
    for si in RSD.s_dictlist:
        data = sokusen_xyzdata(si, RSD)
        x2 += data[0]
        y2 += data[1]
        z2 += data[2]
    z2 = list(np.array(z2)/1000.)
    return(x2,y2,z2)
    #print(sokusen_xyzdata(s_dictlist[0]))

class Sijijiban:
    def __init__(self, input, if4=False, funk1="linear", bai=1, dkey=True, title="Title"):
        s = self
        s.title = title
        #s.my_title = my_title
        #s.file_name = my_title
        s.my_xlabel = "x"
        s.my_ylabel = "y"
        s.input_file = input
        #s.post_file = post
        s.fs = 14
        s.if4 = if4
        s.funk0 = "linear"
        s.funk1 = funk1
        s.my_levels = 0
        s.xb = 0
        s.yb = 0
        s.zb = 0

        s.x = 0
        s.y = 0
        s.z = 0
        s.x_min = 0
        s.x_max = 0
        s.y_min = 0
        s.y_max = 0
        s.z_min = 0
        s.z_max = 0
        s.xi = 0
        s.yi = 0
        s.zi = 0
        s.i_max = 0
        s.i_min = 0
        s.n = 0
        s.bai = bai #データプロットの倍率（4x10=40(x) 40x40=1600プロット
        s.of = 10000 #グラフの描画オフセット
        s.of1 = 100+s.of #範囲外データの作成用オフセット
        if dkey:
            s.set_data() #初期値の更新
        else:
            s.set_data2()
        s.if4_funk() #Trueの場合4点を追加する関数
        s.make_grid() #gridを作成して3次元プロットをする前段

    def set_data(self):
        s = self
        input_file = s.input_file
        # テキストファイルをpandasデータフレーム形式で読み込む
        # 区切り文字はsepで指定できる。例えば、タブ区切りの場合はsep='\t'と記述する
        df = pd.read_csv(input_file, sep=',', encoding="utf-8")
        # x, y, zデータを1次元のnumpyアレイ型へ変換
        x = df['x'].values
        y = df['y'].values
        z = df['z'].values
        s.x = x
        s.y = y
        s.z = z
        s.xb = x
        s.yb = y
        s.zb = z
        #N = len(df) #データ数と最小値、最大値を抽出
        N = 100 #データ数が増えた場合（１万とか）時間がかかりすぎたのでとりあえず100とおいた。
        s.x_min = x.min() # xの最小値
        s.x_max = x.max() # xの最大値
        s.y_min = y.min() # yの最小値
        s.y_max = y.max() # yの最大値
        z_min = z.min()
        z_max = z.max()
        if abs(z_min) > 1000:
            of2 = 2000 #2000mmの拡張
        else:
            of2 = 2 #2mの拡張
        s.i_min = int(z_min-of2) #Zの最小値オフセット array用に整数化
        s.i_max = int(z_max+of2) #Zの最大値
        s.n = N*s.bai #グリッドの数を決定
        
        # 軸の範囲
        s.buf_x = (s.x_max - s.x_min) * 0.001 #元作成者のオフセット
        s.buf_y = (s.y_max - s.y_min) * 0.001 #元作成者のオフセット

    #set data2について、データ追加用関数を用意した法がいいと考えた。後日修正予定220902
    def set_data2(self):
        s = self
        input_file = s.input_file
        # テキストファイルをpandasデータフレーム形式で読み込む
        # 区切り文字はsepで指定できる。例えば、タブ区切りの場合はsep='\t'と記述する
        df = pd.read_csv(input_file, sep=',', encoding="utf-8")
        # x, y, zデータを1次元のnumpyアレイ型へ変換
        x = df['x'].values
        y = df['y'].values
        z = df['z'].values
        da = sokusen_xyz_all()
        #print(type(x),type(da[0]))
        s.xb = np.array(x)
        s.yb = np.array(x)
        s.zb = np.array(x)
        s.x = np.array(list(x) + da[0])
        s.y = np.array(list(y) + da[1])
        s.z = np.array(list(z) + da[2])
        #N = len(df) #データ数と最小値、最大値を抽出
        N = 100 #データ数が増えた場合（１万とか）時間がかかりすぎたのでとりあえず100とおいた。
        s.x_min = x.min() # xの最小値
        s.x_max = x.max() # xの最大値
        s.y_min = y.min() # yの最小値
        s.y_max = y.max() # yの最大値
        z_min = z.min()
        z_max = z.max()
        if abs(z_min) > 1000:
            of2 = 2000 #2000mmの拡張
        else:
            of2 = 2 #2mの拡張
        s.i_min = int(z_min-of2) #Zの最小値オフセット array用に整数化
        s.i_max = int(z_max+of2) #Zの最大値
        s.n = N*s.bai #グリッドの数を決定
        
        # 軸の範囲
        s.buf_x = (s.x_max - s.x_min) * 0.001 #元作成者のオフセット
        s.buf_y = (s.y_max - s.y_min) * 0.001 #元作成者のオフセット

    def if4_funk(self):
        s = self
        if s.if4: #範囲外の4点を追加するパターン(実は最も外の枠のデータをもっと追加したほうが良い可能性あり)
            s.xi = np.array([s.x_min-s.of1,s.x_max+s.of1,s.x_min-s.of1,s.x_max+s.of1]) #追加したい座標
            s.yi = np.array([s.y_min-s.of1,s.y_min-s.of1,s.y_max+s.of1,s.y_max+s.of1]) #追加したい座標
            s.zi = s.makezi(s.x,s.y,s.z,s.xi,s.yi) 
            s.x = ArrayAppend(s.x,s.xi)
            s.y = ArrayAppend(s.y,s.yi)
            s.z = ArrayAppend(s.z,s.zi)

    def make_grid(self):
        s = self
        xi, yi = np.linspace(s.x_min-s.of, s.x_max+s.of, s.n), np.linspace(s.y_min-s.of, s.y_max+s.of, s.n) #グリッド作成用線形補完データ
        xi, yi = np.meshgrid(xi, yi) #グリッドデータの作成
        zi = s.makezi2(s.x,s.y,s.z,xi,yi) #グリッドデータのz線形補完作成
        s.xi = xi
        s.yi = yi
        s.zi = zi

    def get_grid(self,xn,yn): #ほしい位置のzデータを入手したい。
        s = self
        #zi = makezi2(s.x,s.y,s.z,xn,yn) #グリッドデータのz線形補完作成
        zn = s.makezi2(s.x,s.y,s.z,xn,yn)
        return zn

    def makezi(self,x,y,z,xi,yi): #範囲外の4点を作成するためだけの関数(範囲外の補間をするため精度低い)
        s = self
        interp = interpolate.Rbf(x, y, z, epsilon=1,function=s.funk0,mode="1-D") #計算式取得
        zi = interp(xi, yi) #Zメッシュ計算
        return zi

    def makezi2(self,x,y,z,xi,yi): #グリッドの全データを作成する関数 #mt = 'cubic' or 'linear'
        s = self
        zi = interpolate.griddata((x, y), z, (xi, yi), method=s.funk1,rescale=False)
        return zi

def plot_conta(Ama, di=5, my_title="Title", fname="fname",show_boo=False): #コンタ描画（クラスから独立）
    s = Ama
    fs = 6
    zmax = int(( s.i_max - (s.i_max%di))+di) #整数diの倍数で切り上げ。
    zmin = int(( s.i_min - (s.i_min%di))) #整数diの倍数で切り捨て。
    
    my_levels = np.arange(zmin, zmax, di) #diは関数内変数コンタ図の刻み[m or mm]
    plt.rcParams['font.size'] = fs # グラフの基本フォントサイズの設定
    fig = plt.figure(figsize=(18, 12)) # 画像のインスタンスを生成    
    # 等高線
    ax1 = fig.add_subplot(111)
    #contour1 =ax1.contour(s.xi,s.yi,s.zi,5,vmin=-1,vmax=1, levels=my_levels,colors=['black'],)
    contour1 =ax1.contour(s.xi,s.yi,s.zi,vmin=-1,vmax=1, levels=my_levels,colors=['black'],)
    # 等高線の線上に数値ラベルを表示
    #contour1.clabel(fmt='%1.1f', fontsize=fs,inline=False)
    contour1.clabel(fmt='%1.1f', fontsize=fs, inline=True, inline_spacing=1)
    
    ax1.set_title(my_title, fontsize=fs)
    ax1.set_xlabel(s.my_xlabel, fontsize=fs)
    ax1.set_ylabel(s.my_ylabel, fontsize=fs)
    ax1.set_xlim(s.x_min - s.buf_x, s.x_max + s.buf_x)
    ax1.set_ylim(s.y_min - s.buf_y, s.y_max + s.buf_y)
    plt.grid() # 格子グリッド
    plt.gca().set_aspect('equal', adjustable='box') # X,Y軸の1目盛りの表示の縦横比を揃える 
    plt.tick_params(labelsize=fs)
    plt.tight_layout()
    if show_boo:
        plt.show()
    else:
        fig.savefig(fname + "_" + now + ".pdf")
    #plt.show()
    plt.close()
    
def plot_3d(Ama, key, fname): #三次元確認＋画像作成
    s = Ama
    fig = plt.figure(dpi=120)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(s.x_min - s.buf_x, s.x_max + s.buf_x)
    ax.set_ylim(s.y_min - s.buf_y, s.y_max + s.buf_y)
    if s.if4:
        ax.plot_surface(s.xi,s.yi,s.zi,cmap="coolwarm") #面で描画（元）
    else:
        ax.scatter(s.xi,s.yi,s.zi,s=0.1) #点で描画kita
    ax.set_box_aspect((1, 1,1))
    plt.savefig(fname + "3d_"+now+".png",dpi=130)
    if key:
        plt.show()
    else:None
    plt.close()

def plot_3d2(INSL, key, data1, bor_boo=""): #三次元確認＋画像作成
    #s = Ama
    #s2 = Ama2
    Ama0 = INSL[0]

    fig = plt.figure(dpi=120)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(Ama0.x_min - Ama0.buf_x, Ama0.x_max + Ama0.buf_x)
    ax.set_ylim(Ama0.y_min - Ama0.buf_y, Ama0.y_max + Ama0.buf_y)
    
    #ax.plot_surface(s.xi,s.yi,s.zi,cmap="coolwarm") #面で描画（元）
    #ax.scatter(s.x,s.y,s.z,s=5,label="boring survey")
    #ax.scatter(s2.xi,s2.yi,s2.zi,s=0.1,label="Jiyama") #点で描画kita
    for Ama in INSL:
        ax.scatter(Ama.xi,Ama.yi,Ama.zi,s=0.1,label=Ama.title) #点で描画kita
    #ax.set_xlim(Ama0.x_min - Ama0.buf_x, Ama0.x_max + Ama0.buf_x)
    #ax.set_ylim(Ama0.y_min - Ama0.buf_y, Ama0.y_max + Ama0.buf_y)   
    ax.scatter(data1[0], data1[1], data1[2], s=5 ,label="Energy(Hokan)") #点で描画kita
    if bor_boo == "":
        None
    else:
        ax.scatter(bor_boo[0], bor_boo[1], bor_boo[2], s=5 ,label="bor")
    ax.set_box_aspect((1, 1, 1))
    #ax.legend()
    if key:
        plt.show()
    else:None
    plt.close()
    #plt.savefig(fname + "3d2_"+now+".png",dpi=130)
    #if key:
    #    plt.show()
    #else:None

def plot_3d_2data_Ani(AmaL, key, fname, data1): #三次元グラフ回転動画作成
    from matplotlib.animation import FuncAnimation
    from mpl_toolkits.mplot3d import Axes3D

    s = AmaL[0]

    fig = plt.figure(dpi=120)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(s.x_min - s.buf_x, s.x_max + s.buf_x)
    ax.set_ylim(s.y_min - s.buf_y, s.y_max + s.buf_y)

    def plot_graph():
        for Ama in Amal:
            ax.scatter(s.xi,s.yi,s.zi,s=0.1,label=Ama.title) #点で描画kita
            #ax.scatter(s.x,s.y,s.z,s=5,label="boring survey")
        #ax.set_xlim(s2.x_min - s2.buf_x, s2.x_max + s2.buf_x)
        #ax.set_ylim(s2.y_min - s2.buf_y, s2.y_max + s2.buf_y)
        #if s2.if4:
        #    ax.plot_surface(s2.xi,s2.yi,s2.zi,cmap="coolwarm") #面で描画（元）
        #else:
        #    ax.scatter(s2.xi,s2.yi,s2.zi,s=0.1,label="Jiyama") #点で描画kita
        
        ax.scatter(data1[0], data1[1], data1[2], s=5 ,label="Energy(Hokan)") #点で描画kita
        ax.set_box_aspect((1, 1, 1))
        ax.legend()

    def plt_graph3d(angle):
        ax.view_init(azim=angle*2)
    frame_n = 180
    
    ani = FuncAnimation(
        fig,
        func=plt_graph3d,
        frames= frame_n,
        init_func=plot_graph,
        interval=150
    )
    ani.save("rolling5.gif", writer="pillow")
    plt.close()

def post_f_getdata(Ama, RPD): #csvファイルに入力されたデータからx,yを取得してその位置のz作成
    #df = pd.read_csv(input_file, sep=',', encoding="utf-8")
    #x = df['x_input'].values
    #y = df['y_input'].values
    #s_key = df['grid'].values
    x = RPD.x
    y = RPD.y
    
    if RPD.key:
        xi, yi = np.meshgrid(x, y) #mesh作成
    else:
        xi = np.array([x])
        yi = np.array([y])
    zi = Ama.get_grid(xi,yi) ##指定位置のz取得 meshのzデータ作成
    return(xi, yi, zi)

class Read_Post_Data():
    def __init__(self, csv_file):
        df = pd.read_csv(csv_file, sep=',', encoding="ms932")
        self.x = df['x_input'].values
        self.y = df['y_input'].values
        self.s_key = df['grid'].values[0]
        self.z_base = df['基準'].values[0]
        #print()
        if self.s_key == "on":
            key = True
        elif self.s_key == "off":
            key = False
        else:
            print("csv grid on or off error -> s_key = False")
            key = True
        self.key = key

class Read_Section_Data():
    def __init__(self, csv_file):
        self.s_dict = {}
        #csv_file = "sec_測線.csv"
        df = pd.read_csv(csv_file, sep=',', encoding="ms932")
        dfc_list = df.columns.values[1:] #A列を削除

        def BoolFromTF(boo): #FalseはY方向に切断する場合に使う
            if boo == "TRUE":
                return True
            elif boo == "FALSE":
                return False
            else:
                print("error bool string error -> FALSE")
                return False

        #print(dfc_list)
        #print(df["測線1"])
        for dfc in dfc_list:
            si = df[dfc]
            self.s_dict[dfc] = [[float(si[0]), float(si[1])], [float(si[2]), float(si[3])], BoolFromTF(si[4])]
        #print(s_dict)
        #s_dictlist = [i for i in s_dict] #キーのみ取出し
        self.s_dictlist = dfc_list #上に記載のやり方は辞書から取出す場合。今回はその必要もない
    
def post_f_getdata_list(Ama, xylist, key): #xylistからその位置のz作成
    x = xylist[0]
    y = xylist[1]
    if key == True:
        xi, yi = np.meshgrid(x, y) #mesh作成
    else:
        xi = np.array([x])
        yi = np.array([y])
    zi = Ama.get_grid(xi,yi) #meshのzデータ作成
    return xi,yi,zi

def plot_sokusen3(datal, name, fs=14, show=False): #測線位置の地盤高さを確認するだけの関数220518
    plt.rcParams['font.size'] = fs # グラフの基本フォントサイズの設定
    #zi = np.array(zi)*1000.
    # 画像のインスタンスを生成
    fig = plt.figure(figsize=(18, 8))    
    # 等高線
    ax1 = fig.add_subplot(111)
    #ax1.set_xlim(min(yi), max(yi))
    #ax1.set_ylim(0, 40)
    ss = 20
    for datai in datal:
        #print(datai)
        ax1.scatter(datai[4], np.array(datai[2])*1000, s=ss, label=datai[3])
    ax1.set_aspect('equal')
    hans, labs = ax1.get_legend_handles_labels()
    ax1.legend(handles=hans, labels=labs)
    ax1.set_title(name) #タイトルを追加
    # 格子グリッド
    plt.grid()
    # X,Y軸の1目盛りの表示の縦横比を揃える
    #plt.gca().set_aspect('equal', adjustable='box') 
    plt.tick_params(labelsize=fs)
    plt.tight_layout()
    fig.savefig(name + ".jpg")
    if show:
        plt.show()
    plt.close()

#def post_f(Ama,input_file,key,p_name): #csv読み込みと対応するzデータの作成
def post_f(AmaL, p_name, RPD, fs=10): #2次元図面にテキストで深さを記入することを想定
    #以下グラフ化操作
    
    plt.rcParams['font.size'] = fs # グラフの基本フォントサイズの設定
    fig = plt.figure(figsize=(18, 12))  
    ax1 = fig.add_subplot(111)
    
    Ama = AmaL[0]
    z_max = (max(Ama.z))
    z_min = (min(Ama.z))

    for k in range(len(AmaL)):
        Ama = AmaL[k]
        t_offset = 500*k
        xi, yi, zi = post_f_getdata(Ama, RPD)
        zi = zi - RPD.z_base
        ax1.scatter(xi,yi)

        rx = np.min(xi)+(np.max(xi)-np.min(xi))/10
        ry = np.min(yi)+(np.max(yi)-np.min(yi))/10
        ax1.text(rx, ry+t_offset, str(Ama.title), color=plot_color_list[k])
        print(Ama.title)
        
        if max(abs(z_max),abs(z_min)) > 5000: #現状単位判断の手法はスマートではない。
            key = True
        else:
            key = False

        for i in range(len(zi)): #単位がmかmmか判断してから書き込み。
            for j in range(len(zi[0])):
                if key:
                    txt = str("%.0f"%zi[i][j]) #単位mm
                else:
                    txt = str("%.1f"%zi[i][j]) #単位m/10
                ax1.text(xi[i][j], yi[i][j] + t_offset, txt, color=plot_color_list[k])
        
    
    #ax1.set_box_aspect((1, 1,1))
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box') 
    plt.tick_params(labelsize=fs)
    plt.tight_layout()
    fig.savefig("./pdf/" + p_name + ".pdf") #許容できるエラー発生 posx and posy
    fig.savefig("./pdf/" + p_name + ".png")
    plt.close()

def my_dirmake(dir_pass): #保存場所の作成
    new_dir_path_recursive = dir_pass
    try:
        os.makedirs(new_dir_path_recursive)
    except FileExistsError:
        pass

def make_z1():
    input_file = './csv/input.csv' # 読み込むcsvファイル（ここではボーリングデータの支持地盤深さ）
    #funk1 = 'cubic' #'cubic or linear'
    #if4 = False #Falseは範囲内のみ計算,Trueは4点予測（予測の精度は低いので2022年4月時点では使わない。
    #bai = 20 #グリッド計算時の分割数倍数（基データ数x bai） 分割数を増やすと精度もよくなるが、3dで確認するときに動作が遅くなる。
    #bai = 1
    Ama = Sijijiban(input_file, if4=False, funk1='cubic', bai=1, dkey=True, title="oldCubic")
    return(Ama)

def make_z2():
    #input_file = 'Jiyama_manual1.csv'
    input_file = './csv/Jiyama_all.csv'
    #funk1 = 'linear' #'cubic or linear'
    #if4 = False #Falseは範囲内のみ計算,Trueは4点予測（予測の精度は低いので2022年4月時点では使わない。
    #bai = 1 #グリッド計算時の分割数倍数（基データ数x bai） 分割数を増やすと精度もよくなるが、3dで確認するときに動作が遅くなる。
    Ama = Sijijiban(input_file, if4=False, funk1='linear', bai=1, dkey=True, title="Jiyama")
    return(Ama)

def make_z3():
    input_file = './csv/input.csv' # 読み込むcsvファイル（ここではボーリングデータの支持地盤深さ）
    #funk1 = 'cubic' #'cubic or linear'
    #if4 = False #Falseは範囲内のみ計算,Trueは4点予測（予測の精度は低いので2022年4月時点では使わない。
    #bai = 20 #グリッド計算時の分割数倍数（基データ数x bai） 分割数を増やすと精度もよくなるが、3dで確認するときに動作が遅くなる。
    #bai = 1
    Ama = Sijijiban(input_file, if4=False, funk1='cubic', bai=1, dkey=False, title="newCubic")
    return(Ama)

def re_sokusen(data): #data にLのリストを追加 Lの始点はx0,y0とする。
    new = []
    xl = data[0][0]
    yl = data[1][0]
    x0 = xl[0]
    y0 = yl[0]
    #print(xl)
    #print(x0,y0)
    for i in range(len(xl)):
        x1 = xl[i]
        y1 = yl[i]
        #print(x1,y1)
        dx = x1 - x0
        dy = y1 - y0
        L = np.sqrt(dx**2 + dy**2)
        new.append(L)
        #print(L)
    data += tuple([new])
    return(data)

def fname_to_name(name_csv):
    #csvファイル（？）からファイル名に付ける名前の算定
    spncsv = name_csv.split(".")[:-1] #.でスプリットしたリスト
    if len(spncsv) == 1: #リストの長さが１のときはそのままname
            name = spncsv[0]
    elif len(spncsv) == 0: #リスト長さが0のときはそのままname
        name = name_csv
    else:
        name = str.join(spncsv) #リスト長さが2以上の場合は文字列を合体させる
    return(name)

def make_post_only(AmaL, name_csv):
    RPD = Read_Post_Data(name_csv) #csvを読み込んだデータのオブジェクト
    #基礎位置のデータ
    #xi, yi, zi = post_f_getdata(Ama, RPD) #True=grid
    #PDFの作成
    (name, extention) = os.path.splitext(name_csv)
    #name = fname_to_name(name_csv) #csvから.csvの部分を除去したい。
    post_f(AmaL, name, RPD)
    #return([xi, yi, zi])

#post_file = ['post_ランプウェイ.csv', 'post_資源物ヤード.csv', 'post_洗車場棟.csv']
#post_file = ['post_ランプウェイ.csv','post_渡り廊下1.csv','post_渡り廊下2.csv','post_渡り廊下3.csv']
#post_file = ['post_スラグストックヤード棟(2)_new.csv']
#post_file = ['post_渡り廊下2_new.csv']
#grid_if = [True,True,False,True] #グリッドにするかどうか

def make_sokusen_data(pdlist, RSD):
    il = range(len(RSD.s_dictlist))
    #il = [2]
    for i in il: 
        #print(RSD.s_dictlist[i],RSD.s_dict[RSD.s_dictlist[i]])
        #xylist = make_sokusen(s_dict[s_dictlist[i]], 100) #(xylist0, 分割個数）
        da0 = sokusen_xyzdata(RSD.s_dictlist[i], RSD) #測線のdxfから作成されたデータの読み取り。RSDは測線の名前のためだけにも必要
        xylist = [da0[0], da0[1]]
        #data = [0 for i in range(len(pdlist))]
        data = []
        for j in range(len(pdlist)):
            Ama = pdlist[j]
            data1 =  post_f_getdata_list(Ama, xylist, False) #補間データ
            data1 += tuple([Ama.title]) #dataの末尾にタイトルの追加
            data1 = re_sokusen(data1) #dataの最後にLのリストを追加
            data.append(data1)
        #plot_sokusen3は入力データをリスト形式に変更。plot_sokusenシリーズを削除
        #モデルからモデルの断面を作成する。([data, data, data], title)
        plot_sokusen3(data, RSD.s_dictlist[i])  

def make_section_data(pdlist, RSD):
    il = range(len(RSD.s_dictlist)) #H軸, G軸, F軸,..,A軸
    #il = [2]
    for i in il: 
        #print(RSD.s_dictlist[i],RSD.s_dict[RSD.s_dictlist[i]])
        xylist = make_sokusen(RSD.s_dict[RSD.s_dictlist[i]], 100) #(xylist0, 分割個数）2点から間のデータを作成する。
        #da0 = sokusen_xyzdata(RSD.s_dictlist[i], RSD) #測線のdxfから作成されたデータの読み取り。RSDは測線の名前のためだけにも必要
        #xylist = [da0[0], da0[1]]
        #da0 = s_dict[RSD.s_dictlist[i]]
        #xylist = [da0[0], da0[1]]
        #data = [0 for i in range(len(pdlist))]
        data = []
        for j in range(len(pdlist)):
            Ama = pdlist[j]
            data1 =  post_f_getdata_list(Ama, xylist, False) #補間データ
            data1 += tuple([Ama.title]) #dataの末尾にタイトルの追加
            data1 = re_sokusen(data1) #dataの最後にLのリストを追加
            data.append(data1)
        #plot_sokusen3は入力データをリスト形式に変更。plot_sokusenシリーズを削除
        #モデルからモデルの断面を作成する。([data, data, data], title)
        plot_sokusen3(data, RSD.s_dictlist[i])  

def main():
    #保存用フォルダの作成
    my_dirmake("./png")
    my_dirmake("./pdf")
    
    #地盤レベルインスタンスの作成
    #読み込むinput.csvをリスト化しておくのがいい。後日修正予定。
    Ama = make_z1() #補間データ
    #Ama2 = make_z2() #地山データ
    Ama3 = make_z3() #補間データ2
    bor = [Ama.x,Ama.y,Ama.z]

    #3d表示で確認
    name_csv = "./csv/post_エネルギー.csv"
    name_csv0 = "./post_エネルギー.csv"
    RPD = Read_Post_Data(name_csv) #csvを読み込んだデータのオブジェクト
    (name, extention) = os.path.splitext(name_csv)
    #name = fname_to_name(name_csv0) #csvから.csvの部分を除去したい。
    data3 = post_f_getdata(Ama, RPD) #xyデータからAmaモデルで予測される支持層レベル
    plot_3d2([Ama3, Ama], True, data3, bor_boo=bor)
    #各支点の位置のレベルを比較する時に使用する。
    #for csvn in post_file: #各csvファイルを開いて読み込み、xy位置に対応するzデータを作成した後、pdfとpngを作成して保存する。
    #    make_post_only([Ama, Ama2, Ama3], csvn)
    
    #コンタ作成
    conta_Ama = Ama3
    plot_conta(conta_Ama, di = 5, my_title = conta_Ama.title, fname = conta_Ama.title, show_boo=True)
    #断面の描画（測線の断面作成・建物の軸で切断断面作成等）
    #注意を要するのはmake_sokusen　と make_sectionが若干違うこと。dxfの読み取りの有無が違う。
    #RSD = Read_Section_Data("sec_測線.csv") #測線の断面のcsvの読み取り
    #make_sokusen_data([Ama, Ama2, Ama3], RSD)
    #RSD1 = Read_Section_Data("sec_エネ棟軸.csv") #測線の断面のcsvの読み取り
    #make_section_data([Ama, Ama2, Ama3], RSD1)
    
if __name__ == "__main__":
    main()
