# coding=gbk
"""
apriorifunc.py
~~~~~~~~~~~~~~~~
使用马尔科夫链挖掘数据
输入：dataarray 行为模式库中的频繁1集
      minsupmap   各频繁1集的最小支持度
输出：
行为模式库中满足最小支持度和最小置信度的强关联规则
可能的提示错误：
minsupmap目前全部指定0.1
最小置信度目前全部指定0.04  之后需要根据实际情况修改
"""
import pymysql
import time
import datetime

ITEMNUM = 30000
CODESNUM=1
SUPPORT = 0.1
CONFIDENCE = 0.1
SEQLEN=4

class Data:
    def __init__(self, _bv):
        self.bv = _bv  # 行为序号
        self.tidlist = [0] * ITEMNUM  # 支持项位向量
        self.sup = 0.0  # 支持度

    def show(self):
        print("bv=%d:tidlist=%s:sup=%f" % (self.bv, self.tidlist, self.sup))


class Code:
    def __init__(self, _id, _code):
        self.id = _id
        self.code = _code


class Stack:
    """模拟栈"""

    def __init__(self):
        self.items = []

    def isEmpty(self):
        return len(self.items) == 0

    def push(self, item):
        self.items.append(item)

    def pop(self):
        return self.items.pop()

    def peek(self):
        if not self.isEmpty():
            return self.items[len(self.items) - 1]

    def size(self):
        return len(self.items)

    def show(self):
        for item in self.items:
            print item,


class Large:
    def __init__(self):
        self.p = [] #保存频繁集各项元素
        self.sup = []   #保存该频繁集的支持度

def getConnection():
    conn = pymysql.connect(host='localhost',db='supervision', user='root', passwd='root', port=3306, charset='utf8')#之后可以放到配置文件中读取
    return conn
#特殊挖掘情况保存频繁集
def saveProcessLar(Lar,target):
    conn=getConnection()
    cur=conn.cursor()
    sql="select count(*) from information_schema.columns where table_schema='supervision' and table_name='data_process_processinfo_middle';"
    cur.execute(sql)
    res = cur.fetchall()
    num=res[0][0]-6#减去id support user,begintime,endtime和createtime
    if len(Lar)>num:
        for j in range (num,len(Lar)):
            sql="alter table "
            sql+=target['middle']
            sql+=" add `process%d` "%(j+1)
            sql+=" varchar(255) DEFAULT NULL"
            cur.execute(sql)
    ###
    basesql="INSERT INTO " + target['middle'] + "("
    codes=target['codes']
    for ii in range(1,len(Lar)):
        Large=Lar[ii]
        for index in range(len(Large.p)):
            sql=basesql
            for ind in range(len(Large.p[index])):
                sql+="`process%d`, "%(ind+1)
            sql+="`support`,`user`,`begin_time`,`end_time`"
            sql+=") VALUES("
            for ind in range(len(Large.p[index])):
                sql+="'%s', "%(codes[Large.p[index][ind]])
            sql+=str(Large.sup[index])
            sql+=", "
            sql+="'%s' ,"%(target['user'])
            sql+="'%s' ,"%(target['begin_time'])
            sql+="'%s'"%(target['end_time'])
            sql+=")"
            cur.execute(sql)
    cur.close()
    conn.commit()
    conn.close()
def saveLar(Lar,target):
    conn=getConnection()
    cur=conn.cursor()
    sql="TRUNCATE TABLE `"
    sql+=target['middle']
    sql+="`;"
    cur.execute(sql)
    basesql="INSERT INTO " + target['middle'] + "("
    basesql2 = "SELECT * FROM "
    basesql2 += target['codes']
    basesql2 += " where id="
    field=target['field']
    # sql += str(row[i])
    # cur.execute(sql)
    # res = cur.fetchall()
    # rule += [(res[0][1])]
    for Large in Lar:
        for index in range(len(Large.p)):
            if isinstance(Large.p[index],int):#针对Lar[0]特殊处理
                sql=basesql2+str(Large.p[index])
                cur.execute(sql)
                res = cur.fetchall()
                sql=basesql
                sql+="`%s` ,`support` ) VALUES("%(field[res[0][2]])
                sql+="'%s'"%(res[0][1])
                sql+=" , "
                sql+=str(Large.sup[index])
                sql+=")"
                cur.execute(sql)
            else:
                rule=[]
                sql=basesql
                for ind in range(len(Large.p[index])):
                    sql2=basesql2+str(Large.p[index][ind])
                    cur.execute(sql2)
                    res = cur.fetchall()
                    sql+="`%s`, "%(field[res[0][2]])
                    rule+=[res[0][1]]
                sql+="`support` "
                sql+=") VALUES("
                for ind in range(len(Large.p[index])):
                    sql+="'%s', "%(rule[ind])
                sql+=str(Large.sup[index])
                sql+=")"
                cur.execute(sql)
    cur.close()
    conn.commit()
    conn.close()


# tidlist与操作
def tidlistandfunc(tidlist1, tidlist2):
    ret = [0] * ITEMNUM
    for i in range(ITEMNUM):
        if 1 == tidlist1[i] and 1 == tidlist2[i]:
            ret[i] = 1
        else:
            ret[i] = 0
    return ret


# tidlist计数有多少个1
def tidlistCountfunc(tidlist):
    count = 0
    for i in tidlist:
        if i != 0:
            count += 1
    return count


# 图
class Graph(object):
    def __init__(self, dataarray, minsupmap, Lar):
        self.node_neighbors = {}
        self.sup = {}
        self.tidlist = {}
        self.add_nodes(dataarray)
        L2 = Large()
        for i in range(len(dataarray) - 1):
            j = i + 1
            while j < len(dataarray):
                vi = dataarray[i]
                vj = dataarray[j]
                sup1 = min(vi.sup, vj.sup)
                sup2 = max(minsupmap[vi.bv], minsupmap[vj.bv])
                if sup1 >= sup2:
                    andtidlist = tidlistandfunc(vi.tidlist, vj.tidlist)
                    count = tidlistCountfunc(andtidlist)
                    sup3 = (float)(count) / ITEMNUM
                    if sup3 >= sup2:
                        L2.p += [(vi.bv, vj.bv)]
                        L2.sup += [sup3]
                        self.add_edge((vi.bv, vj.bv))
                j += 1
        Lar += [L2]

    def getfirstnode(self):
        return sorted(self.node_neighbors.items(), key=lambda n: (-(self.sup[n[0]]),n[0]))[0][0]

    def add_nodes(self, dataarray):
        for node in range(len(dataarray)):
            self.add_node(dataarray[node].bv, dataarray[node].sup, dataarray[node].tidlist)

    def add_node(self, node, _sup, _tidlist):
        if not node in self.nodes():
            self.node_neighbors[node] = []
            self.sup[node] = _sup
            self.tidlist[node] = _tidlist

    def add_edge(self, edge):
        u, v = edge
        if (v not in self.node_neighbors[u]):
            self.node_neighbors[u].append(v)

    def nodes(self):
        return self.node_neighbors.keys()

    def deletenode(self):
        a = self.getfirstnode()
        self.node_neighbors.pop(a)

    def displaygraph(self):
        arcs = sorted(self.node_neighbors.items(), key=lambda n: len(n[1]), reverse=True)
        for arc in arcs:
            print arc[0]
            print "~~~"
            print arc[1]
            print "**********"
            # 深度优先搜索来获取频繁n集

    def dfssearch(self, minsupmap, Lar):
        s = Stack()
        ps = Stack()
        p = self.getfirstnode()
        s.push(p)
        i = 0
        p = self.node_neighbors[s.peek()][i] if len(self.node_neighbors[s.peek()]) > 0 else -1
        while not s.isEmpty():
            while p > 0:
                if s.size() >= 2:
                    allbeside = True
                    for sf in s.items:
                        beside = False
                        for pf in self.node_neighbors[sf]:
                            if pf == p:
                                beside = True
                                break
                        if not beside:
                            allbeside = False
                            break
                    if allbeside:
                        sup1 = self.sup[p]
                        sup2 = minsupmap[p]
                        for item in s.items:
                            sup1 = min(sup1, self.sup[item])
                            sup2 = max(sup2, minsupmap[item])
                        if sup1 >= sup2:
                            bit = self.tidlist[p]
                            for item in s.items:
                                bit = tidlistandfunc(bit, self.tidlist[item])
                            sup3 = (float)(tidlistCountfunc(bit)) / ITEMNUM
                            if sup3 >= sup2:
                                if len(Lar) == s.size() - 1 + 1:
                                    Lar += [Large()]
                                nodes = []
                                nodes += s.items
                                nodes += [p]
                                Lar[s.size()].p += [tuple(nodes)]
                                Lar[s.size()].sup += [sup3]
                                s.push(p)
                                ps.push(i)
                                i = 0
                            else:
                                i += 1
                        else:
                            i += 1
                    else:
                        i += 1
                else:
                    s.push(p)
                    ps.push(i)
                    i = 0
                p = self.node_neighbors[s.peek()][i] if len(self.node_neighbors[s.peek()]) > i else -1
            else:
                s.pop()
                if not ps.isEmpty():
                    i = ps.pop()
                    i += 1
                    p = self.node_neighbors[s.peek()][i] if len(self.node_neighbors[s.peek()]) > i else -1
###频繁集中的规律添加至数据库
# def autoInsert(rules_table, frontfield, backfield, constantfield, frontdata, backdata, constantdata):
#     sql = "INSERT IGNORE INTO " + rules_table + "("
#     for i in range(len(constantfield)):
#         sql += "`" + constantfield[i] + "`,"
#     for i in range(len(frontdata)):
#         if frontdata[i]>0:
#             sql += "`" + frontfield[i] + "`,"
#     for i in range(len(backdata)):
#         if backdata[i]>0:
#             sql += "`" + backfield[i] + "`,"
#     sql = sql[:-1]
#     sql += ") VALUES("
#     for i in range(len(constantdata)):
#         sql += "%lf," % (constantdata[i])
#     for i in range(len(frontdata)):
#         if frontdata[i]>0:
#             sql += "%d," % (frontdata[i])
#     for i in range(len(backdata)):
#         if backdata[i]>0:
#             sql += "%d," % (backdata[i])
#     sql = sql[:-1]
#     sql += ")"
#     return sql
###挖掘频繁集中的规律
# def rules(Lar, target):
#     conn =getConnection()
#     cur = conn.cursor()
#     for i in range(2, len(Lar) + 1):
#         for j in range(len(Lar[i - 1].sup)):
#             n = 1
#             for h in range(i):
#                 n *= 2
#             for h in range(1, n - 1):
#                 temp1 = []
#                 temp2 = []
#                 tp1 = 0
#                 tp2 = 0
#                 t = 1
#                 counts = 0.0
#                 for m in range(i):
#                     if (h / t) % 2 == 1:
#                         temp1 += [Lar[i - 1].p[j][m]]
#                         tp1 += 1
#                     else:
#                         temp2 += [Lar[i - 1].p[j][m]]
#                         tp2 += 1
#                     t *= 2
#                 for jj in range(len(Lar[tp1 - 1].sup)):
#                     eq = True
#                     for k in range(tp1):
#                         if tp1 == 1:
#                             if Lar[tp1 - 1].p[jj] != temp1[0]:
#                                 eq = False
#                                 break
#                         else:
#                             if Lar[tp1 - 1].p[jj][k] != temp1[k]:
#                                 eq = False
#                                 break
#                     if eq:
#                         counts = Lar[tp1 - 1].sup[jj] * ITEMNUM
#                         break
#                 if Lar[i - 1].sup[j] * ITEMNUM / counts >= CONFIDENCE:  # 检查是否大于等于最小置信度阈值 目前指定0.3
#                     seq1=[-1]*SEQLEN
#                     seq2=[-1]*SEQLEN
#                     basesql="SELECT `column` FROM "
#                     basesql+=target['count']
#                     basesql+=" where `pid`="
#                     for idx in range(len(temp1)):
#                         sql=basesql
#                         sql+=str(temp1[idx])
#                         cur.execute(sql)
#                         res=cur.fetchall()
#                         seq1[res[0][0]]=temp1[idx]
#                     for idx in range(len(temp2)):
#                         sql=basesql
#                         sql+=str(temp2[idx])
#                         cur.execute(sql)
#                         res=cur.fetchall()
#                         seq2[res[0][0]]=temp2[idx]
#                     print temp1,
#                     print "==>",
#                     print temp2
#                     print "置信度:",
#                     confidence = (Lar[i - 1].sup[j] * ITEMNUM / counts)
#                     print confidence * 100,
#                     print "%"
#                     print "支持度",
#                     support = (Lar[i - 1].sup[j])
#                     print support * 100,
#                     print "%"
#                     constantdata = [support, confidence]
#                     sql = autoInsert(target['rules_table'], target['frontfield'], target['backfield'],
#                                      target['constantfield'], seq1, seq2, constantdata)
#                     cur.execute(sql)
#     cur.close()
#     conn.commit()
#     conn.close()

# 获得数据
def get_dataarray(target, targetcount, field):
    dataarray = []
    conn =getConnection()
    cur = conn.cursor()
    sql= "SELECT COUNT(*) FROM "
    sql+= targetcount
    cur.execute(sql)
    res = cur.fetchall()
    global CODESNUM
    CODESNUM=res[0][0]
    sql = "SELECT pid FROM "
    sql += targetcount
    sql += " where (sum/%d)>%f order by sum desc,`pid`" % (ITEMNUM, SUPPORT)
    cur.execute(sql)
    res = cur.fetchall()
    for row in res:
        d = Data(row[0])
        dataarray += [d]
    for d in dataarray:
        count = 0.0
        sql = "SELECT id FROM "
        sql += target
        sql += " where"
        for i in range(len(field)):
            sql += "`%s`" % (field[i])
            sql += "=%d or" % (d.bv)
        sql = sql[:-3]
        cur.execute(sql)
        res = cur.fetchall()
        for row in res:
            d.tidlist[row[0] - 1] = 1
            count += 1
        d.sup = count / ITEMNUM
    cur.close()
    conn.commit()
    conn.close()
    return dataarray
def get_ITEMNUM(target):
    conn =getConnection()
    cur = conn.cursor()
    sql = "SELECT * FROM "
    sql += target
    cur.execute(sql)
    res = cur.fetchall()
    global ITEMNUM
    ITEMNUM = len(res)
    cur.close()
    conn.commit()
    conn.close()
# 获得强关联规则
def getrules(dataarray, minsupmap, target):
    Lar = []
    L1 = Large()
    for d in dataarray:
        L1.p += [d.bv]
        L1.sup += [d.sup]
    Lar += [L1]
    g = Graph(dataarray, minsupmap, Lar)
    while len(g.nodes()) > 0:
        g.dfssearch(minsupmap, Lar)
        g.deletenode()
    if 1==target['type']:
        saveLar(Lar,target)
        #rules(Lar, target)#挖掘关联规则 目前不需要
    else:
        saveProcessLar(Lar,target)
#调试用 显示规则中文名
# def translate(target):
#     conn =getConnection()
#     cur = conn.cursor()
#     sql = "SELECT * FROM "
#     sql += target['rules_table']
#     cur.execute(sql)
#     res = cur.fetchall()
#     for row in res:
#         rule = []
#         column=[]
#         for i in range(1, len(row) - 3):
#             if row[i]:
#                 sql = "SELECT * FROM "
#                 sql += target['codes']
#                 sql += " where id="
#                 sql += str(row[i])
#                 cur.execute(sql)
#                 res = cur.fetchall()
#                 rule += [(res[0][1])]
#                 column+=[(res[0][2])]
#         # rule += [str(row[-2])]
#         # rule += [str(row[-1])]
#         print str(rule).decode('string_escape')
#         sql = "INSERT IGNORE INTO " + target['result'] + "("
#         for i in range(len(column)):
#             sql+="`%s` ,"%(target['field'][column[i]])
#         sql+="`support` "
#         sql += ") VALUES("
#         for i in range(len(rule)):
#             sql += "'%s'," % (rule[i])
#         sql += "'%s'" % (str(row[-3]))
#         sql += ")"
#         cur.execute(sql)
#     cur.close()
#     conn.commit()
#     conn.close()
# 获得规则
def get_table_rules(name, field):
    target = {}
    frontfield = []
    backfield = []
    constantfield = ["support", "confidence"]
    global SEQLEN
    SEQLEN=len(field)
    for i in range(0, len(field)):
        frontfield += ["front_" + field[i]]
        backfield += ["back_" + field[i]]
    #rules_table = name + "_apriorirules"
    middle=name+"_middle"
    target['field']=field
    target['frontfield'] = frontfield
    target['backfield'] = backfield
    #target['rules_table'] = rules_table
    target['middle']=middle
    target['constantfield'] = constantfield
    get_ITEMNUM(name)
    dataarray = get_dataarray(name, name + "_count", field)
    target['codes'] = name + "_codes"
    target['count']=name+"_count"
    target['type']=1
    minsupmap = {}
    i = CODESNUM
    while (i > 0):
        minsupmap[i] = SUPPORT
        i -= 1
    getrules(dataarray, minsupmap, target)
    #translate(target)  #调试用 显示规则中文名
###数据清洗函数###
def filter(middletable,ruletable,field,selectfield,num):#num为取得条数，num为-1时取所有条
    conn=getConnection()
    cur=conn.cursor()
    sql="SELECT "
    for f in selectfield:
        sql+="`%s`, "%(f)
    sql+="`support`,`create_time` FROM "
    sql+=middletable
    sql+=" where "
    for ff in field:
        if ff in selectfield:
            sql+=" `%s` is not null and"%(ff)
        else:
            sql+=" `%s` is null and"%(ff)
    sql=sql[:-3]
    sql+=" order by `support` desc"
    if num>0:
        sql+=" limit %d"%(num)
    cur.execute(sql)
    res = cur.fetchall()
    for r in res:
        sql="Insert into %s ("%(ruletable)
        for sf in selectfield:
            sql+="`%s` ,"%(sf)
        sql+="`support`,`create_time`,`type`) VALUES("
        for i in range(len(selectfield)+2):#包括support和create_time
            sql+="'%s',"%(r[i])
        sql+="%d)"%(len(selectfield))
        cur.execute(sql)
    cur.close()
    conn.commit()
    conn.close()
def get_ip_packet_rules():
    field=["time","host_id","send_mac_address","recv_mac_address","send_ip","send_port","recv_ip","recv_port"]
    get_table_rules("data_process_ippacket", field)
    middle="data_process_ippacket_middle"
    result="data_process_ippacket_result"
    ###频繁三项集###
    filter(middle,result,field,["time", "host_id", "recv_ip"],-1)
    ###频繁二项集###
    filter(middle,result,field,[ "host_id", "time"],5)
    filter(middle,result,field,[ "host_id", "recv_ip"],5)
    filter(middle,result,field,[ "time", "recv_ip"],5)
def get_warning_information_rules():
    field = ["time", "userid",  "description","rank", "species"]
    get_table_rules("warning_information", field)
    # middle="warning_information_middle"
    # result="warning_information_result"
    # ###频繁三项集###
    # filter(middle,result,field,["time", "userid", "description"],-1)
    # ###频繁二项集###
    # filter(middle,result,field,[ "userid", "description"],5)
    # filter(middle,result,field,[ "userid", "rank"],5)
    # filter(middle,result,field,[ "userid", "species"],5)
def get_data_process_fileinfo_file_rules():
    field = ["time", "file_name", "user", "operate_type", "host_id"]
    get_table_rules("data_process_fileinfo_file", field)
    middle="data_process_fileinfo_file_middle"
    result="data_process_fileinfo_file_result"
    ###频繁四项集###
    filter(middle,result,field,["time", "file_name", "user", "operate_type"],-1)
    ###频繁三项集###
    filter(middle,result,field,["time", "file_name", "user"],5)
    filter(middle,result,field,["time",  "user", "operate_type"],5)
    filter(middle,result,field,["time", "file_name", "operate_type"],5)
    filter(middle,result,field,["file_name", "user", "operate_type"],5)
    ###频繁二项集###
    filter(middle,result,field,["user", "operate_type"],10)
    filter(middle,result,field,["file_name", "user"],10)
def get_data_process_fileinfo_type_rules():
    field = ["time", "file_type", "user", "operate_type", "host_id"]
    get_table_rules("data_process_fileinfo_type", field)
    middle="data_process_fileinfo_type_middle"
    result="data_process_fileinfo_type_result"
    ###频繁四项集###
    filter(middle,result,field,["time", "file_type", "user", "operate_type"],-1)
    ###频繁三项集###
    filter(middle,result,field,["time", "file_type", "user"],5)
    filter(middle,result,field,["time",  "user", "operate_type"],5)
    filter(middle,result,field,["time", "file_type", "operate_type"],5)
    filter(middle,result,field,["file_type", "user", "operate_type"],5)
    ###频繁二项集###
    filter(middle,result,field,["user", "operate_type"],10)
    filter(middle,result,field,["file_type", "user"],10)
def get_data_process_resource_warning_rules():
    field = ["time", "user", "process_name", "resources_name", "warning_rank"]
    get_table_rules("data_process_resource_warning", field)
    # middle="data_process_resource_warning_middle"
    # result="data_process_resource_warning_result"
    # ###频繁五项集###
    # filter(middle,result,field,["time", "user", "process_name", "resource_name", "warning_rank"],-1)
    # ###频繁四项集###
    # filter(middle,result,field,["time", "user", "process_name", "resource_name"],5)
    # filter(middle,result,field,["time", "user", "process_name", "warning_rank"],5)
    # filter(middle,result,field,["time", "user", "resource_name", "warning_rank"],5)
    # filter(middle,result,field,["time", "process_name", "resource_name", "warning_rank"],5)
    # filter(middle,result,field,["user", "process_name", "resource_name", "warning_rank"],5)
    # ###频繁三项集###
    # filter(middle,result,field,[ "user", "process_name", "resource_name"],5)
    # filter(middle,result,field,[ "user", "process_name", "warning_rank"],5)
    # filter(middle,result,field,["time", "user", "process_name"],5)
    # filter(middle,result,field,["time", "process_name", "resource_name"],5)
    # filter(middle,result,field,["user", "process_name", "warning_rank"],5)
    # filter(middle,result,field,["time", "user", "warning_rank"],5)
    # ###频繁二项集###
    # filter(middle,result,field,[ "user", "process_name"],5)
    # filter(middle,result,field,[ "process_name", "resource_name"],5)
    # filter(middle,result,field,[ "process_name", "warning_rank"],5)
    # filter(middle,result,field,["time","process_name"],5)
    # filter(middle,result,field,["time","resource_name"],5)
    # filter(middle,result,field,["user","resource_name"],5)
#特殊挖掘
def get_data_process_processinfo_rules(date,begin_time,end_time):
    global ITEMNUM
    ITEMNUM = 4
    global SUPPORT
    SUPPORT = 0.75
    begin_dates=[]
    end_dates=[]
    begin_date=datetime.datetime.strptime(date+" "+begin_time, '%Y-%m-%d %H:%M:%S')
    begin_dates+=[begin_date]
    end_date=datetime.datetime.strptime(date+" "+end_time, '%Y-%m-%d %H:%M:%S')
    end_dates+=[end_date]
    for d in range(ITEMNUM-1):
        delta = datetime.timedelta(days=1)
        begin_date=begin_date+delta
        begin_dates+=[begin_date]
        end_date=end_date+delta
        end_dates+=[end_date]
    conn =getConnection()
    cur = conn.cursor()
    sql="DROP TABLE IF EXISTS `"
    sql+="data_process_processinfo_middle"
    sql+="`;"
    cur.execute(sql)
    ###
    sql="CREATE TABLE if not exists "
    sql+="data_process_processinfo_middle"
    sql+=" (`id` int(11) NOT NULL AUTO_INCREMENT,`create_time` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,`support` double NOT NULL ,"
    sql+= "`user` varchar(45) NOT NULL,`begin_time` varchar(45) NOT NULL,`end_time`  varchar(45) NOT NULL , "
    sql+=" PRIMARY KEY (`id`)) ENGINE=InnoDB DEFAULT CHARSET=utf8;"
    cur.execute(sql)
    ###
    users=[]
    sql="SELECT DISTINCT `user` FROM supervision.data_process_processinfo_raw;"
    cur.execute(sql)
    res = cur.fetchall()
    for row in res:
        users+=[row[0]]
    for user in users:
        processlist = [0] * ITEMNUM
        processlist[0] = []
        processlist[1] = []
        processlist[2] = []
        processlist[3] = []
        for dd in range(ITEMNUM):
            sql = "SELECT `process_name` FROM supervision.data_process_processinfo_raw where `begintime`>='%s' and `endtime`<='%s' and `user`='%s';"%(begin_dates[dd],end_dates[dd],user)
            cur.execute(sql)
            res = cur.fetchall()
            for row in res:
                processlist[dd] += [row[0]]
            processlist[dd] = list(set(processlist[dd]))
        Codesmap = {}
        for plist in processlist:
            for p in plist:
                Codesmap[p] = 1
        Codes = []
        id = 1
        for key in Codesmap.keys():
            code = Code(id, key)
            Codes += [code]
            id += 1
        global CODESNUM
        CODESNUM=len(Codes)
        Codesmap.clear()
        for c in Codes:
            Codesmap[c.id] = c.code
        dataarray = []
        for c in Codes:
            sum = 0.0
            data = Data(c.id)
            for i in range(len(processlist)):
                for j in range(len(processlist[i])):
                    if c.code == processlist[i][j]:
                        sum += 1
                        data.tidlist[i] = 1
                        processlist[i][j] = c.id
            if sum / ITEMNUM >= SUPPORT:
                data.sup = sum / ITEMNUM
                dataarray += [data]
        dataarray = sorted(dataarray,key=lambda d:(-(d.sup),d.bv))
        minsupmap = {}
        i = CODESNUM
        while (i > 0):
            minsupmap[i] = SUPPORT
            i -= 1
        target={}
        target['type']=2
        target['codes']=Codesmap
        target['middle']="data_process_processinfo_middle"
        target['user']=user
        target['begin_time']=begin_time
        target['end_time']=end_time
        getrules(dataarray,minsupmap,target)
    sql="DROP TABLE IF EXISTS data_process_processinfo_result;"
    cur.execute(sql)
    sql="DROP TABLE IF EXISTS data_process_processinfo_result_a;"
    cur.execute(sql)
    sql="DROP TABLE IF EXISTS data_process_processinfo_result_b;"
    cur.execute(sql)
    sql+="CREATE TABLE data_process_processinfo_result_a LIKE data_process_processinfo_middle ; "
    sql+="CREATE TABLE data_process_processinfo_result_b LIKE data_process_processinfo_middle  ;"
    sql+="CREATE TABLE data_process_processinfo_result LIKE data_process_processinfo_middle  ;"
    sql+="ALTER TABLE data_process_processinfo_result ADD `type` int(11);"
    cur.execute(sql)
    sql="select count(*) from information_schema.columns where table_schema='supervision' and table_name='data_process_processinfo_middle';"
    cur.execute(sql)
    res = cur.fetchall()
    num=res[0][0]-6#减去id support user,begintime,endtime和createtime
    ### a ###
    sql="insert INTO data_process_processinfo_result_a SELECT * FROM supervision.data_process_processinfo_middle order by "
    for n in range(num,0,-1):
        sql+=" `process%d` desc,"%(n)
    sql=sql[:-1]
    sql+=" limit 10"
    cur.execute(sql)
    sql="ALTER TABLE data_process_processinfo_result_a ADD `type` int(11) DEFAULT 1;"
    cur.execute(sql)
    ###从middle表中取出a中的元素# ###
    sql="SELECT id FROM data_process_processinfo_result_a"
    cur.execute(sql)
    res=cur.fetchall()
    for r in res:
        sql="delete from data_process_processinfo_middle where id= %d"%(r[0])
        cur.execute(sql)
    ### b ###
    sql="insert INTO data_process_processinfo_result_b SELECT * FROM supervision.data_process_processinfo_middle order by `support` desc,"
    for n in range(num,0,-1):
        sql+=" `process%d` desc,"%(n)
    sql=sql[:-1]
    sql+=" limit 15"
    cur.execute(sql)
    sql="ALTER TABLE data_process_processinfo_result_b ADD `type` int(11) DEFAULT 2;"
    cur.execute(sql)
    sql="Insert into data_process_processinfo_result select * from data_process_processinfo_result_a"
    cur.execute(sql)
    sql="Insert into data_process_processinfo_result select * from data_process_processinfo_result_b"
    cur.execute(sql)
    sql="DROP TABLE data_process_processinfo_result_a"
    cur.execute(sql)
    sql="DROP TABLE data_process_processinfo_result_b"
    cur.execute(sql)
    cur.close()
    conn.commit()
    conn.close()


if __name__ == '__main__':
    start = time.clock()
    #get_ip_packet_rules()
    get_warning_information_rules()
    #get_data_process_fileinfo_file_rules()
    #get_data_process_fileinfo_type_rules()
    #get_data_process_resource_warning_rules()
    #get_data_process_processinfo_rules('2017-03-03','10:00:00','11:00:00')
    end = time.clock()
    print str(end - start) + "s"
