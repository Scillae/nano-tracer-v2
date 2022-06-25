from infra.Reader import Reader

def strands_construct(data, top_path, traj_path):
    varname = 'strands'
    # load data from result_catch
    SL_result_catch(data, varname, 'load')
    if data['SL_content']:
        catched_result = data['SL_content']
        data['SL_content'] = None
        return catched_result
    
    
    # not catched
    reader = Reader(top_path, traj_path)
    strands_series = reader.read_data()


    # save in result_catch
    data['SL_content'] = strands_series
    SL_result_catch(data, varname, 'save')
    return strands_series