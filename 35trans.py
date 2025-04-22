#-----------------------------------------------------------------------------------------
#----------------------------区域内为可微调列表--------------------------------------------------
#核酸序列 
l_base=['u','a','c','g','t','n'] #lowercase
u_base=['U','A','C','T','G',  #下方为IUPAC所命名指代
        'N','M','B','Y','H','R','S','K','D','W']#uppercase
n_residue=l_base+u_base#全部碱基

#程序默认修饰基团在aa右边，如有基团在左，请添加在下方l表格中
#核酸修饰基团（左）
l=['V','f','25r','o','x','m','p','5m','xd','sg','rg','fx','C12','C6','SB','Cy3#','r','C3']
#核酸修饰集团（link）
la=['*','#','.','•','^','-','Ps','s']
#核酸修饰基团（右），精准切片
r=['MOE','DCA','Chol','idAB','iSp9','3ddC','TriGalNAc']
residue=n_residue + l + r + la 

#strange(没有意义的标点符号)
dev=['[',']','/',"'",'(',')',' ','′',"'"]
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
def section_1 (Seq_orig):
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

#0.3 序列 终·标准化【将所有非aa基团剔除:得到 “ Seq_norm ”】
def ult_uni (Seq_uni):
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
    for i in range(0,len(Seq_list)):
    #若末位为aa残基，删除
        if Seq_list[-1] in n_residue:
            Seq_list=Seq_list[:-1]
        else:
    #若末位非aa残基，终·标准化，并统计残基个数(即position)
            Seq_u=ult_uni(uni_seq(Seq_list))
            location=len(Seq_u)
        #这里判断修饰左,linkage或是右
            if Seq_list[-1] in l:
                l_seq=("{}:{}".format(Seq_list[-1],location))
                modi.append(l_seq)
            elif Seq_list[-1] in la:    
                la_seq=("{}:{}->{}".format(Seq_list[-1],location,location+1))
                modi.append(la_seq)               
            else:
                m_seq=("{}:{}".format(Seq_list[-1],location+1))
                modi.append(m_seq)
            Seq_list=Seq_list[:-1]
    #有时候我们需要区分lowercase/uppercase
    #那我们重新跑一次,输出所有碱基位置
    for i in range(0,len(seq_list)):
        if seq_list[-1] in n_residue:
            seq_u=ult_uni(uni_seq(seq_list))
            u_seq=("{}:{}".format((seq_list[-1]),len(seq_u)))
            modi.append(u_seq)
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
    Seq_norm=ult_uni(uni_seq(clear(section_1(Seq_orig))))
    print(f'Seq_norm:{''.join(Seq_norm)}')
#4.1 输出修饰表格    
    Seq_list=clear(section_1(Seq_orig))
    con_table(con_dict(seek(Seq_list)))
#----------------------------------这是一条装订线-------------------------------------------------
#----------------------------------------------------------------------------------------------------




