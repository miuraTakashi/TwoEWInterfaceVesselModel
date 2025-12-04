"""
TwoEWInterfaceVesselModel - 血管壁の動きのシミュレーションモデル

血管壁の動きをシミュレーションする数値計算モデルです。
2つの血管の界面の時間発展を、拡散、ランダムノイズ、弾性相互作用を考慮して計算します。
"""

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from copy import deepcopy
from matplotlib import animation
import os
from pathlib import Path
import datetime


# 血管壁の動きの支配方程式：
# 
# ∂_t h_1 = d_h ∂_x^2 h_1 + σ η_1(x,t) - k_s/w_1(h_1-(H-w_1)) + k_e/(2r)(h_2-h_1+2r)
# ∂_t h_2 = d_h ∂_x^2 h_2 + σ η_2(x,t) - k_e/(2r)(h_2-h_1+2r) + k_s/w_2(h_3-h_2+w_2)
# ∂_t h_3 = d_h ∂_x^2 h_3 + σ η_3(x,t) - k_s/w_2(h_3-h_2+w_2) + k_e/(2r)(h_4-h_3+2r)
# ∂_t h_4 = d_h ∂_x^2 h_4 + σ η_4(x,t) - k_e/(2r)(h_4-h_3+2r) - k_s/w_1(h_4+H-w_1)
#
# パラメータの説明:
# - h_i(x,t): 時刻tで場所xのi番目の血管の界面の位置
# - d_h: 血管壁の表面張力。炎症時には細胞間接着が外れて減少する
# - σ: 内皮細胞のランダムな運動の強さ
# - η_i(x,t): 平均0、分散1の正規分布の乱数（標準正規乱数）
# - k_s: 周辺の結合組織の硬さ（単位長さあたりのバネ係数）
# - k_e: 血管壁同士の相互作用の単位長さあたりのバネ係数
# - H: 血管系全体の高さの範囲
# - w_1: 血管と表皮の間の間隔
# - w_2: 血管同士の間の間隔
# - r: 血管の半径


def TwoVesselModel(L=5, dh=1, sigma=1, ks=0.1, ke=0.1, w1=1, w2=2, r=0.5, fps=20, T=1000, progress=True):
    """
    2つの血管の界面位置の時間発展を計算する関数
    
    パラメータ:
    -----------
    L : float
        空間領域の長さ（デフォルト: 5）
    dh : float
        表面張力係数 d_h（デフォルト: 1）
    sigma : float
        ノイズ強度 σ（デフォルト: 1）
    ks : float
        結合組織の硬さ k_s（デフォルト: 0.1）
    ke : float
        血管間相互作用 k_e（デフォルト: 0.1）
    w1 : float
        血管と表皮の間隔（デフォルト: 1）
    w2 : float
        血管間の間隔（デフォルト: 2）
    r : float
        血管の半径（デフォルト: 0.5）
    fps : int
        アニメーションのフレームレート（デフォルト: 20）
    T : float
        シミュレーション時間（デフォルト: 1000）
    progress : bool
        進捗バーの表示（デフォルト: True）
    
    戻り値:
    -------
    hList : list
        各時刻での界面位置 [h1, h2, h3, h4] のリスト
    """
    H = w2/2 + 2*r + w1
    print(H)
    dx = 0.1
    dt = 0.001
    sqrt_dt = np.sqrt(dt)

    n = int(L/dx)
    loop = int(T/dt)
    loopPerOutput = int(T/dt/100)

    def Delta(h):
        """2階微分を計算（周期境界条件）"""
        return (np.roll(h, 1) + np.roll(h, -1) - 2*h) / dx / dx

    hList = []
    h1 = np.zeros(n) + w2/2 + 2*r
    h1New = deepcopy(h1)
    h2 = np.zeros(n) + w2/2
    h2New = deepcopy(h2)
    h3 = np.zeros(n) - w2/2
    h3New = deepcopy(h3)
    h4 = np.zeros(n) - w2/2 - 2*r
    h4New = deepcopy(h4)

    iterator = tqdm(range(loop)) if progress else range(loop)
    
    for t in iterator:
        h1New += (dt * (
            -ks/w1*(h1-(H-w1)) + ke/(2*r)*(h2-h1+2*r)
            + dh*Delta(h1)
        ) + sigma * np.random.randn(n) * sqrt_dt)

        h2New += (dt * (
            -ke/(2*r)*(h2-h1+2*r) + ks/w2*(h3-h2+w2)
            + dh*Delta(h2)
        ) + sigma * np.random.randn(n) * sqrt_dt)

        h3New += (dt * (
            -ks/w2*(h3-h2+w2) + ke/(2*r)*(h4-h3+2*r)
            + dh*Delta(h3)
        ) + sigma * np.random.randn(n) * sqrt_dt)

        h4New += (dt * (
            -ke/(2*r)*(h4-h3+2*r) - ks/w1*(h4+H-w1)
            + dh*Delta(h4)
        ) + sigma * np.random.randn(n) * sqrt_dt)

        h1 = deepcopy(h1New)
        h1[0] = h1[-1] = w2/2 + 2*r
        h1[h1 > H] = H
        h2 = deepcopy(h2New)
        h2[0] = h2[-1] = w2/2
        h3 = deepcopy(h3New)
        h3[0] = h3[-1] = -w2/2
        h4 = deepcopy(h4New)
        h4[h4 < -H] = -H
        h4[0] = h4[-1] = -w2/2 - 2*r

        if t % loopPerOutput == 0:
            hList.append([h1, h2, h3, h4])
    
    return hList


def plotVessel(h, ax=None, save_path=None):
    """
    血管の界面位置をプロットする関数
    
    パラメータ:
    -----------
    h : list
        [h1, h2, h3, h4] のリスト
    ax : matplotlib.axes.Axes, optional
        プロットする軸。指定しない場合は新しい図を作成
    save_path : str, optional
        保存先のパス。指定しない場合はresultsフォルダに自動保存
    """
    if ax is None:
        fig, ax = plt.subplots()
        create_new_figure = True
    else:
        create_new_figure = False
    
    ax.plot(h[0], color='red')
    ax.plot(h[1], color='red')
    ax.plot(h[2], color='red')
    ax.plot(h[3], color='red')
    ax.set_ylim(-5, 5)
    ax.set_aspect('equal')
    ax.fill_between(range(len(h[0])), h[0], h[1], color='red', alpha=0.2)
    ax.fill_between(range(len(h[2])), h[2], h[3], color='red', alpha=0.2)
    
    if create_new_figure:
        # resultsフォルダを作成
        results_dir = Path("results")
        results_dir.mkdir(exist_ok=True)
        
        # 保存パスを決定
        if save_path is None:
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            save_path = results_dir / f"vessel_{timestamp}.png"
        else:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
        
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"画像を保存しました: {save_path}")


def export_vessel_animation(dh=0.1, sigma=1, ks=0.001, ke=0.1, w1=1, w2=2, r=0.5, T=100, fps=20):
    """
    引数で与えた dh, sigma, ks, ke を用いてシミュレーションし、
    形態変化のアニメーションを MP4 で保存する。
    保存ファイル名はパラメータ入りで自動命名。
    
    パラメータ:
    -----------
    dh : float
        表面張力係数（デフォルト: 0.1）
    sigma : float
        ノイズ強度（デフォルト: 1）
    ks : float
        結合組織の硬さ（デフォルト: 0.001）
    ke : float
        血管間相互作用（デフォルト: 0.1）
    w1 : float
        血管と表皮の間隔（デフォルト: 1）
    w2 : float
        血管間の間隔（デフォルト: 2）
    r : float
        血管の半径（デフォルト: 0.5）
    T : float
        シミュレーション時間（デフォルト: 100）
    fps : int
        フレームレート（デフォルト: 20）
    
    戻り値:
    -------
    filename : str
        保存されたファイル名
    """
    # インタラクティブ表示を抑止
    prev_interactive = plt.isinteractive()
    plt.ioff()
    try:
        # 進捗バーも抑止
        result = TwoVesselModel(dh=dh, sigma=sigma, ks=ks, ke=ke, w1=w1, w2=w2, r=r, T=T, fps=20, progress=False)

        # resultsフォルダを作成
        results_dir = Path("results")
        results_dir.mkdir(exist_ok=True)
        
        filename = results_dir / f"vessel_dh{dh}_sigma{sigma}_ks{ks}_ke{ke}.mp4"
        fig, ax = plt.subplots()

        def animate(i):
            ax.clear()
            plotVessel(result[i], ax=ax)
            ax.set_title(f"Frame {i}")

        ani = animation.FuncAnimation(fig, animate, frames=len(result), interval=1000.0/fps, repeat=False)

        try:
            from matplotlib.animation import FFMpegWriter
            writer = FFMpegWriter(fps=fps, metadata=dict(artist=''), bitrate=1800)
            ani.save(filename, writer=writer)
        except Exception:
            ani.save(filename, writer='ffmpeg', fps=fps)
        finally:
            plt.close(fig)
    finally:
        plt.interactive(prev_interactive)

    return filename


if __name__ == "__main__":
    # 使用例1: 基本的なシミュレーション
    print("=== 基本的なシミュレーション ===")
    result = TwoVesselModel(L=5, dh=0.1, sigma=1, ks=0.1, ke=0.1, T=100)
    plotVessel(result[-1], save_path="results/vessel_basic.png")
    
    # 使用例2: 異なるパラメータでの比較
    print("\n=== 異なるパラメータでの比較 ===")
    for ks in [10, 0]:
        result = TwoVesselModel(ks=ks, ke=1, sigma=1, T=100)
        print(f"ks={ks}")
        plotVessel(result[-1], save_path=f"results/vessel_ks{ks}.png")
    
    # 使用例3: アニメーションの生成
    print("\n=== アニメーションの生成 ===")
    filename = export_vessel_animation(dh=0.1, sigma=1, ks=0.1, ke=0.1)
    print(f"アニメーションを保存しました: {filename}")
    
    # 使用例4: 統計解析（コメントアウト - 時間がかかるため）
    # print("\n=== 統計解析 ===")
    # resultList = []
    # for i in range(100):
    #     result = TwoVesselModel(T=100)
    #     resultList.append(result)
    # 
    # h1hatList = []
    # h2hatList = []
    # h3hatList = []
    # h4hatList = []
    # for i in range(100):
    #     h1hatList.append(np.abs(np.fft.fft(resultList[i][-1][0])))
    #     h2hatList.append(np.abs(np.fft.fft(resultList[i][-1][1])))
    #     h3hatList.append(np.abs(np.fft.fft(resultList[i][-1][2])))
    #     h4hatList.append(np.abs(np.fft.fft(resultList[i][-1][3])))
    # 
    # h1hatMean = np.array(h1hatList).mean(axis=0)
    # h2hatMean = np.array(h2hatList).mean(axis=0)
    # h3hatMean = np.array(h3hatList).mean(axis=0)
    # h4hatMean = np.array(h4hatList).mean(axis=0)
    # 
    # fig, ax = plt.subplots()
    # ax.plot(h1hatMean)
    # ax.plot(h2hatMean)
    # ax.plot(h3hatMean)
    # ax.plot(h4hatMean)
    # ax.set_xlim(1, 25)
    # ax.set_ylim(0.1, 10)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.legend(["h1", "h2", "h3", "h4"])
    # 
    # results_dir = Path("results")
    # results_dir.mkdir(exist_ok=True)
    # plt.savefig(results_dir / "fourier_analysis.png", dpi=150, bbox_inches='tight')
    # plt.close(fig)
    # print(f"統計解析結果を保存しました: {results_dir / 'fourier_analysis.png'}")

