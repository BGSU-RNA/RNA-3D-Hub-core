from functions import *


def main():
    """
    """
    (cursor, conn) = connect_to_mysql()    
    select_database(cursor)    
    create_motifs_table(cursor)
    create_history_table(cursor)    
    
    loops, groups = read_motifs_csv_file()
    m = transform_data(loops, groups)
    newVersion, oldVersion = get_version(cursor)

    if oldVersion is not None:
        L, G = get_previous_release(cursor,oldVersion)    
        M = transform_data(L, G)
        newMotifs, history = compare_releases(cursor, m, M, loops, L)
        load_motifs_data(cursor, conn, newMotifs, newVersion)
        load_history_data(cursor, conn, history, newVersion)
    else:
        pdb.set_trace()
        for group,loops in m.iteritems():
            newId = generate_id(cursor,loops)
            del m[group]
            m[newId] = loops
        load_motifs_data(cursor, conn, m, newVersion)
    
#     pdb.set_trace()
                
        
    close_mysql_connection(cursor,conn)    
     
    print 'OK'


if __name__ == '__main__':
    main()
    
    
    
    
#         load_data_array(cursor,conn,loops,groups,newVersion)     