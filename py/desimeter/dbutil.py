

def dbquery(comm,operation,parameters=None) :

    cx=comm.cursor()
    cx.execute(operation,parameters)
    names=[d.name for d in cx.description]
    res=cx.fetchall()
    cx.close()

    results=dict()
    for i,name in enumerate(names) :
        results[name]=[r[i] for r in res]

    return results

def get_petal_ids(comm) :

    cx=comm.cursor()
    cx.execute("select relname from pg_class where relkind='r' and relname !~ '^(pg_|sql_)';")
    tables=[d[0] for d in cx.fetchall()]
    cx.close()

    # I am sure one can do better
    petalids=[]
    for table in tables :
        if table.find("positioner_moves_p")>=0 :
            tmp=table.replace("positioner_moves_p","")
            #print(tmp)
            try :
                i=int(tmp)
                if i<20 :
                    print(table)
                    petalids.append(i)
            except ValueError :
                pass
    return petalids

def get_pos_ids(comm, petal_id) :
    res=dbquery(comm,"select distinct pos_id from posmovedb.positioner_moves_p%s",(int(petal_id),))
    return res["pos_id"]

def get_petal_loc(petal_id) :
    # hardcoded locations of petal_id, should go to data
    petal_id2loc = dict()
    petal_id2loc[2]=7
    petal_id2loc[3]=3
    petal_id2loc[4]=0
    petal_id2loc[5]=1
    petal_id2loc[6]=2
    petal_id2loc[7]=8
    petal_id2loc[8]=4
    petal_id2loc[9]=9
    petal_id2loc[10]=5
    petal_id2loc[11]=6
    # lbl petal 1 tests
    petal_id2loc[1]=0

    return petal_id2loc[petal_id]
