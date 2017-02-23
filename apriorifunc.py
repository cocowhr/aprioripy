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
ITEMNUM和ITEMLEN目前需要手动指定 可改成自动获取数据库行为模式库大小和行为模式数据长度
minsupmap目前全部指定0.1
最小置信度目前全部指定0.04  之后需要根据实际情况修改
"""
import pymysql

ITEMNUM = 1285


class Data:
    def __init__(self, _bv):
        self.bv = _bv  # 行为序号
        self.tidlist = [0] * ITEMNUM  # tidlist
        self.sup = 0.0  # 支持度

    def show(self):
        print("bv=%d:tidlist=%s:sup=%f" % (self.bv, self.tidlist, self.sup))


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
        self.p = []
        self.sup = []


def printlar(Lar):
    for Large in Lar:
        for index in range(len(Large.p)):
            print Large.p[index], Large.sup[index]


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
        return sorted(self.node_neighbors.items(), key=lambda n: len(n[1]), reverse=True)[0][0]

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
        p = self.node_neighbors[p][i] if len(self.node_neighbors[p]) > 0 else -1
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
                        p = self.node_neighbors[p][i] if len(self.node_neighbors[p]) > 0 else -1
                    else:
                        i += 1
                else:
                    s.push(p)
                    ps.push(i)
                    i = 0
                    p = self.node_neighbors[p][i] if len(self.node_neighbors[p]) > 0 else -1
            else:
                s.pop()
                if not ps.isEmpty():
                    i = ps.pop()
                    i += 1
                    p = self.node_neighbors[s.peek()][i] if len(self.node_neighbors[s.peek()]) > i else -1


def rules(Lar):
    for i in range(2, len(Lar) + 1):
        j = 0
        for j in range(len(Lar[i - 1].sup)):
            n = 1
            h = 0
            for h in range(i):
                n *= 2
            h = 0
            for h in range(1, n - 1):
                temp1 = []
                temp2 = []
                tp1 = 0
                tp2 = 0
                t = 1
                counts = 0.0
                m = 0
                for m in range(i):
                    if (h / t) % 2 == 1:
                        temp1 += [Lar[i - 1].p[j][m]]
                        tp1 += 1
                    else:
                        temp2 += [Lar[i - 1].p[j][m]]
                        tp2 += 1
                    t *= 2
                jj = 0
                for jj in range(len(Lar[tp1 - 1].sup)):
                    eq = True
                    k = 0
                    for k in range(tp1):
                        if tp1 == 1:
                            if Lar[tp1 - 1].p[jj] != temp1[0]:
                                eq = False
                                break
                        else:
                            if Lar[tp1 - 1].p[jj][k] != temp1[k]:
                                eq = False
                                break
                    if eq:
                        counts = Lar[tp1 - 1].sup[jj] * ITEMNUM
                        break
                if Lar[i - 1].sup[j] * ITEMNUM / counts >= 0.04:  # 检查是否大于等于最小置信度阈值 目前指定0.04
                    print temp1,
                    print "==>",
                    print temp2
                    print "置信度:",
                    print (Lar[i - 1].sup[j] * ITEMNUM / counts) * 100,
                    print "%"
                    print "支持度",
                    print (Lar[i - 1].sup[j]) * 100,
                    print "%"


# 获得数据
def getdataarray():
    dataarray = []
    try:
        conn = pymysql.connect(host='localhost', user='root', passwd='root', port=3306, charset='utf8')
        cur = conn.cursor()
        cur.execute("USE apriori")
        cur.execute("SELECT * FROM apriori.seq")
        res=cur.fetchall()
        ITEMNUM=len(res)
        cur.execute("SELECT pid FROM apriori.count where sum/(select count(*) from apriori.seq)>0.1 order by sum desc;")
        res = cur.fetchall()
        for row in res:
            d = Data(row[0])
            dataarray += [d]
        for d in dataarray:
            count = 0.0
            cur.execute("SELECT id FROM apriori.seq where protocol=%d or process=%d or url=%d or time=%d;" % (
                d.bv, d.bv, d.bv, d.bv))
            res = cur.fetchall()
            for row in res:
                d.tidlist[row[0] - 1] = 1
                count += 1
            d.sup = count / ITEMNUM
        cur.close()
        conn.commit()
        conn.close()
    except  Exception:
        print("error")
    return dataarray

#获得强关联规则
def getrules(dataarray, minsupmap):
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
    rules(Lar)


if __name__ == '__main__':
    dataarray = getdataarray()
    minsupmap = {}
    i = 50
    while (i > 0):
        minsupmap[i] = 0.1
        i -= 1
    getrules(dataarray,minsupmap)