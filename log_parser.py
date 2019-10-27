import pandas as pd

def find_lines(lines,key,variable_type,delim=':'):
    parsed=None
    for line in lines:
        if key == line.split(delim)[0]:
            label=line.split(delim)[0]
            data=line.split(delim)[1]
            if variable_type==0:
                data=[data.strip()]
            elif variable_type==1:
                data=data.split()
                data=[i.strip() for i in data]
            else:
                print("varialbe_type of",key,"is wrong")
            parsed=(label,data)
            break
        elif key=='Ratio' and 'Ratio < 0' in line[:10]:
            parsed=('Ratio',['0'])
    return parsed

def parse_file(path,key_vtype_tuple_list):
    f=open(path,'r')
    f_lines=f.readlines()
    f.close()
    parsed=[]
    for key_vtype_tuple in key_vtype_tuple_list:
        key=key_vtype_tuple[0]
        vtype=key_vtype_tuple[1]
        find_result=find_lines(f_lines,key,vtype)
        if find_result==None:
            raise NameError('Error in',key,path)
        parsed.append(find_result)
    return parsed

def file_todf(log_files,key_vtype_tuple_list,dtype):
    cname=None
    table=[]
    for log_file in log_files:
        cname_temp=[]
        category=None
        try:
            key_data_tuple_list=parse_file(log_file,key_vtype_tuple_list)
        except Exception as e:
            print(log_file,"passed",e)
            continue
        row=[]
        for (key,data) in key_data_tuple_list: # see if there is cateogory data
            if 'cate' in key or 'Cate' in key:
                category=data
        for (key,data) in key_data_tuple_list: # parsing a file to cname and data
            row=row+data
            if len(data)==1: #key of single value
                cname_temp=cname_temp+[key]
            else: #key of multiple values
                #print(len(category),key,len(data))
                if category==None:
                    for i in range(len(data)):
                        cname_temp=cname_temp+[key+'.'+str(i)]
                elif len(category)==len(data):
                    for i in range(len(data)):
                        cname_temp=cname_temp+[key+'.'+category[i]]
                else:#maybe something have gone wrong
                    print("category found but not mathing data")
        if cname==None or cname==cname_temp:
            cname=cname_temp
        else:
            print("column not matching")
        table.append([log_file]+row)
    cname=['filename']+cname
    return pd.DataFrame(table,columns=cname,dtype=dtype)