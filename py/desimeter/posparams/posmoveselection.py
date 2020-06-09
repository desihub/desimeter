import numpy as np

def posmove_selection(table,log_note_selection=None) :
    """
    select entries in input table and return table with selected rows
    """
    if log_note_selection is None : return table
    print("Unique values of LOG_NOTE before selection: \n{}".format(np.unique(table["LOG_NOTE"])))
    selection_strs = log_note_selection.split('|')
    ok = np.array([False for i in range(len(table))])
    for s in selection_strs:
        ok |= get_matches(table, s)
    print("Unique values of LOG_NOTE after selection: \n{}".format(np.unique(table["LOG_NOTE"][ok])))
    print("Number of selected entries = {}".format(np.sum(ok)))

    return table[ok]

def get_matches(table, selection_str):
    keywords=[]
    for v in selection_str.split("&") :
        v = v.strip()
        if len(v)>0 : keywords.append(v)
    print("Posmove selection keywords: {}".format(keywords))
    ok=np.repeat(True,len(table))
    for keyword in keywords :
        for i in range(len(table)) :
            ok[i] &= (str(table["LOG_NOTE"][i]).find(keyword)>=0)
    return ok