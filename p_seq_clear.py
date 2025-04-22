#-----------------------------------------------------------------------------------------
#----------------------------区域内为可微调列表--------------------------------------------------
#氨基酸残基 
aa_s_upper=['Q','W','A','R','N','C','E','G','H','D',
      'I','L','K','M','F','P','S','T','Y','V','X']#Single
aa_s_lower=['r','q','f','k','v','p']#有的文献中用lowercase表示D-amino_acid
aa_s_substitute=['O','Φ','Ω','φ','U','Nal','3-Quin','(K/Q)','4gp','KLH','SKB',
                 'Xaa1','Xaa2','Xaa3']#Biopython的seq1函数很笨，不能读低于三字母的转换,这里是单字母手动转
aa_substitute=['N-Me-bAla','a-MeLys','aMeTyr','W(7-Me)','aMeGlu','Asp(OMe)','F(4-2ae)','K(Ac)',
               'α-MeLeu','hTyr','(2-Ant)Ala','Lys(iPr)','(3-I)Tyr','Lys(Ac)','Lys(CysAcid-DOTA(Ga))',
               'Phe[4-(2-aminoethoxy)]','a-MePhe','Phe(4-OMe)','Lys(COCF3)','Lys(COtBu)','W(7-Ph)','Orn(COMe)']#需视为整体的带修饰氨基酸（不转X）
aa_s=aa_s_lower+aa_s_upper+aa_s_substitute
aa_t=['Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile',
      'Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val']#Ternate
aa_u=['Nle','Cit','2-Nal','Hse','2Nal','1Nal','Abu','Pen','Cha',
      'Ahx','Hyp','AzK','Aib','Pgy','Dml','Pia','cit','THP','Dab','Orn','Achc']#uncommon
aa_trans=aa_t+aa_u#需转换残基
aa_residue=aa_s+aa_t+aa_u#Seq_norm切片参照物
aa=aa_residue+aa_substitute

#程序默认修饰基团在aa右边，如有基团在左，请添加在下方l表格中
#修饰基团(左)
l=['5FAM','(D)','o','m','7-Ph','a-','a-Me','4-(2-aminoethoxy)','Me','formyl','cyclo','aMe',
   'N-Me-b','d','D-','(3-I)','(isoindole)','(2-Ant)','Ac','Propionic_acid','Propionic_acid',
   '3,3,3-Trifluoropropionic acid','Pentanoic acid','pr','N3_Acid','FPrpTriazoleMe_Acid']
#修饰基团(linkage)
la=['s','linker']
#修饰基团(右)[精准切片]
r=['N3','OH','*','NH2','PEG2','CPQ2','Me','SKB','KLH','Ac','hex-5-ynoyl','PEG4','PEG12','(7-Me)',
   '(4-2ae)','Q-tag','(OMe)','CH2F','(iPr)','(Ac)','Acid','(Ga-BL34)','(CysAcid-DOTA(Ga))','(CysAcid-amido-N,N-dimethyl-ammoniomethyl-trifluoroborate)',
   '(CysAcid-triazole-N,N-dimethyl-ammoniomethyl-trifluoroborate)','(PepBF3-BL40)','(AMBF3-BL41)','4-OMe',
   'COCF3','COtBu','ac','COMe']
residue=aa_residue + l + r + la
section_2=residue+aa_substitute

#strange(没有意义的标点符号)
dev=['[',']','/',"'",'(',')',' ','′',"'",'-',' ']
#---------------------------------------------------------------------------------------
#------------------该区域为调用的第三方库，有些需pip install-------------------------------
#使用seq1函数将aa简写转换
from Bio.SeqUtils import seq1
#制表
from itertools import zip_longest
#关闭程序前给您问个好
import time
#----------------------------------------------------------------------------------------
#---------------------------该区域为定义函数合集-------------------------------------------
#0.0 序列粗切片【切出各种已知残基：氨基酸残基 & 已知基团】
#    (definition：粗标准化序列（Seq_list）为可操作序列)
def section_1 (Seq_orig,residue):
    Seq_list=[]
    i=0
    n=len(Seq_orig)
    while i < n:
        matched = False
        #优先检查更长的匹配项
        for length in range(min(len(max(residue, key=len)), n - i), 0, -1):
            current = Seq_orig[i:i+length]
            if current in residue:
                Seq_list.append(current)
                i += length
                matched = True
                break
        if not matched:
            Seq_list.append(Seq_orig[i])
            i += 1
    return Seq_list

#0.1 序列粗清除（清除各种无意义标点）
def clear(Seq_list):
    for i in range(len(Seq_list)-1, -1, -1):
        if Seq_list[i] in dev:
            del Seq_list[i]
    return Seq_list

#0.2 序列粗标准化【氨基酸残基转换：缩写和转X】
def uni_seq(Seq_list):
    Seq_uni=[]
    for i in Seq_list:
        if i in aa_trans:
            Seq_uni.append(seq1(i))
        elif i in aa_s_lower:
            Seq_uni.append(i.upper())
        elif i in aa_s_substitute:
            Seq_uni.append('X')
        else:    
            Seq_uni.append(i)
    return Seq_uni

#0.3 序列 终·标准化【将所有非aa基团剔除:得到 “ Seq_norm ”】
def ult_uni (Seq_uni,aa_s_upper):
    Seq_norm=[]
    Seq_norm=[i for i in Seq_uni if i in set(aa_s_upper)]
    return Seq_norm

#1.0 找到修饰基团
def seek (Seq_list,aa):
#1.1 首先建立一个空列表装修饰位点结果
    modi=[]
#1.2 从末位开始遍历切片列表(Seq_uni)
    #不止遍历一次，我们给另一次遍历赋值
    seq_list=Seq_list
    for i in range(0,len(Seq_list)):
    #若末位为aa残基，删除
        if Seq_list[-1] in aa:
            Seq_list=Seq_list[:-1]
        else:
    #若末位非aa残基，终·标准化，并统计残基个数(即position)
            Seq_u=ult_uni(uni_seq(Seq_list),aa)
            location=len(Seq_u)
        #这里判断修饰左,linkage或是右
            if Seq_list[-1] in l:
                l_seq=("{}:{}".format(Seq_list[-1],location+1))
                modi.append(l_seq)
            elif Seq_list[-1] in la:    
                la_seq=("{}:{}->{}".format(Seq_list[-1],location,location+1))
                modi.append(la_seq)               
            else:
                m_seq=("{}:{}".format(Seq_list[-1],location))
                modi.append(m_seq)
            Seq_list=Seq_list[:-1]
    #以上忽略了没有角标的X和非常见氨基酸
    #那我们为了这些X重新跑一次
    for i in range(0,len(seq_list)):
        #定位原文为‘X’的氨基酸
        if seq_list[-1]=='X':
            seq_u=ult_uni(uni_seq(seq_list),aa)
            X_seq=("{}:{}".format((seq_list[-1]),len(seq_u)))
            modi.append(X_seq)
            seq_list=seq_list[:-1]
        #定位原文为uncommon的氨基酸
        elif seq_list[-1] in aa_u:
            seq_u=ult_uni(uni_seq(seq_list),aa)
            u_seq=("{}:{}".format((seq_list[-1]),len(seq_u)))
            modi.append(u_seq)
            seq_list=seq_list[:-1]
        #定位原文为小写指代的氨基酸
        elif seq_list[-1] in aa_s_lower:
            seq_u=ult_uni(uni_seq(seq_list),aa)
            X_seq=("{}:{}".format((seq_list[-1]),len(seq_u)))
            modi.append(X_seq)
            seq_list=seq_list[:-1]
        #定位原文已转‘X’的氨基酸
        elif seq_list[-1] in aa_s_substitute:
            seq_u=ult_uni(uni_seq(seq_list),aa)
            X_seq=("{}:{}".format((seq_list[-1]),len(seq_u)))
            modi.append(X_seq)
            seq_list=seq_list[:-1]        
        #定位原文带修饰的且需视作整体的氨基酸
        elif seq_list[-1] in aa_substitute:
            seq_u=ult_uni(uni_seq(seq_list),aa)
            m_seq=("{}:{}".format((seq_list[-1]),len(seq_u)))
            modi.append(m_seq)
            seq_list=seq_list[:-1]
        else:
            seq_list=seq_list[:-1]
    return modi

#2.0 建立修饰字典
def con_dict(modi):
    modi_dict={}
    for i in modi:
        key=i.split(':',1)[0]
        if key in modi_dict:
            modi_dict[key].append(i)
        else:
            modi_dict[key]=[i]
    return modi_dict

#3.0 将字典内容表格化，方便复制粘贴
def con_table(modi_dict):
    data={key: [i.split(':')[1] for i in value] for key, value in modi_dict.items()}
    # 填充缺失值为 None
    max_len = max(len(v) for v in data.values())
    filled_data = {k: list(v) + [] * (max_len - len(v)) for k, v in data.items()}
    print(f'已检测到修饰有：{len(filled_data)}种')
    for key, values in filled_data.items():
        row = [key] + values
        print(",".join(f"{str(item):^1}" for item in row))
    return
#-------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------程序在此运行--------------------------------------------------------------------------------------------
while True:
    Seq_orig=input("请粘贴完整蛋白序列(输入1可退出):")
    if Seq_orig=='1':
        print('辛苦啦~')
        time.sleep(2.5)
        break
#4.0 输出Seq_norm
    Seq_norm=ult_uni(uni_seq(clear(section_1(Seq_orig,residue))),aa_s_upper)
    print(f'Seq_norm:{''.join(Seq_norm)}')
#4.1 输出修饰表格    
    Seq_list=clear(section_1(Seq_orig,section_2))
    con_table(con_dict(seek(Seq_list,aa)))

#----------------------------------这是一条装订线-------------------------------------------------
#----------------------------------------------------------------------------------------------------







