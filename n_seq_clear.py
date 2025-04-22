#-----------------------------------------------------------------------------------------
#----------------------------区域内为可微调列表--------------------------------------------------
#核酸序列 
l_base=['u','a','c','g','t'] #lowercase
u_base=['U','A','C','T','G',  #下方为IUPAC所命名指代
        'N','B','Y','H','R','S','K','D','W','V','β','δ']#uppercase
substitute_base=['oC','oU','mC','(invdT)','(invdA)','Ghd','Ahd','Chd','Uhd','Tam']#需要替换的特殊修饰碱基（不转n）
n_residue=l_base+u_base#Seq_norm切片参照物
base=n_residue+substitute_base#全部碱基

#程序默认修饰基团在base右边，如有基团在左，请添加在下方l列表中
#核酸修饰基团（左）
l=['P','25r','m','TexasRed','a-Me','Cy3','d','(Cy3)#','o','x','xd','(NAG25)','(NAG25)s',
   '(NAG31)s','p','5m','rg','fx','C12','C6','SB','Cy3#','r','FAM','Yak','m5','m6','(NAG37)',
   '(NAG31)s(invAb)s','(NAG25)s(invAb)s','(NAG37)s(invAb)s','(NAG37)s','(NAG25)s(invAb)']
#核酸修饰集团（link）
la=['*','#','.','•','^','-','Ps','s','p*']
#核酸修饰基团（右），精准切片
r=['MOE','DCA','Chol','idAB','iSp9','3ddC','TriGalNAc','C3','BHQ2','BHQ1','IBFQ','(SB)','Agn'
   ,'DCA','f','e','(uPACE)','(cPACE)','L96','M','(invAb)','s(invAb)','hd']
residue=n_residue + l + r + la
section_2=residue+substitute_base 

#strange(没有意义的标点符号)
dev=['[',']','/',"'",'(',')',' ','′',"'",' ']
#---------------------------------------------------------------------------------------
#------------------该区域为调用的第三方库，有些需pip install-------------------------------
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

#0.2 序列粗标准化【碱基残基转换：小写和转n】
def uni_seq(Seq_list):
    Seq_uni=[]
    for i in Seq_list:
        if i in n_residue:
            i=i.lower()
            if i not in l_base:
                Seq_uni.append('n')
            else :
                Seq_uni.append(i)
        else:
            Seq_uni.append(i)
    return Seq_uni

#0.3 序列 终·标准化【将所有非l_base基团剔除:得到 “ Seq_norm ”】
def ult_uni (Seq_uni,l_base):
    Seq_norm=[]
    Seq_norm=[i for i in Seq_uni if i in set(l_base)]
    return Seq_norm

#1.0 找到修饰基团
def seek (Seq_list):
#1.1 首先建立一个空列表装修饰位点结果
    modi=[]
#1.2 从末位开始遍历切片列表(Seq_uni)
    #不止遍历一次，我们给另一次遍历赋值
    seq_list=Seq_list
    SEQ_list=Seq_list
    for i in range(0,len(Seq_list)):
    #若末位为碱基，删除
        if Seq_list[-1] in base:
            Seq_list=Seq_list[:-1]
        else:
    #若末位碱基，剔除所有非碱基元素，并统计残基个数(即position)
            Seq_u=ult_uni(Seq_list,base)
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
    #有时候我们需要区分lowercase/uppercase/或特殊含碱基修饰指代
    #那我们重新跑一次,输出所有需要的碱基位置
    for i in range(0,len(seq_list)):
        if seq_list[-1] in base:
            seq_u=ult_uni(seq_list,base)
            u_seq=("{}:{}".format((seq_list[-1]),len(seq_u)))
            modi.append(u_seq)
            seq_list=seq_list[:-1]
        else:
            seq_list=seq_list[:-1]
    return modi

#1.3 直接给出lowercase/uppercase,方便复制
def Up_or_Low (Seq_orig):
    modi=[]
    SEQ_list=ult_uni(clear(section_1 (Seq_orig,residue)),n_residue) 
    for i in range(0,len(SEQ_list)):
        if SEQ_list[-1].isupper():
            seq_upper=('{}:{}'.format('uppercase',len(SEQ_list)))
            modi.append(seq_upper)
            SEQ_list=SEQ_list[:-1]
        elif SEQ_list[-1].islower():
            seq_lower=('{}:{}'.format('lowercase',len(SEQ_list)))
            modi.append(seq_lower)
            SEQ_list=SEQ_list[:-1]
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
    for key, values in filled_data.items():
        row = [key] + values
        print(",".join(f"{str(item):^1}" for item in row))
    return
#-------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------程序在此运行--------------------------------------------------------------------------------------------
while True:
    Seq_orig=input("请粘贴完整核酸序列(输入1可退出):")
    if Seq_orig=='1':
        print('辛苦啦~')
        time.sleep(2.5)
        break
#4.0 输出Seq_norm
    Seq_norm=ult_uni(uni_seq(clear(section_1(Seq_orig,residue))),l_base)
    print(f'Seq_norm:{''.join(Seq_norm)}')
#4.1 输出修饰表格    
    Seq_list=clear(section_1(Seq_orig,section_2))
    con_table(con_dict(seek(Seq_list)))
    con_table(con_dict(Up_or_Low(Seq_orig)))

#----------------------------------这是一条装订线-------------------------------------------------
#----------------------------------------------------------------------------------------------------







