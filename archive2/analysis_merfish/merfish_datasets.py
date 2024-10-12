merfish_datasets = {
    # exp 1-1
    'P14NRa_ant': "merfish_06142023/ant/region0",
    'P28NRa_ant': "merfish_06142023/ant/region1",
    
    # exp 1-2
    'P14NRa_pos': "merfish_06142023/pos/region0",
    'P28NRa_pos': "merfish_06142023/pos/region1",

    # exp 2-1
    'P21NRa_ant': "merfish_20231114/region0",
    'P21DRa_ant': "merfish_20231114/region2",
    'P28DRa_ant': "merfish_20231114/region1",
    
    # exp 2-2
    'P21NRa_pos': "merfish_20231120/region0",
    'P21DRa_pos': "merfish_20231120/region1",
    'P28DRa_pos': "merfish_20231120/region2",
    
    # exp 3-1
    'P14NRb_ant': "merfish_202404051211/region_0", 
    'P28NRb_ant': "merfish_202404051211/region_1",
    
    # exp 3-2
    'P14NRb_pos': "merfish_202404051214/region_0", 
    'P28NRb_pos': "merfish_202404051214/region_1",
    
    # exp 4-1
    'P21NRb_ant': "merfish_202404091526/region_2",
    'P21DRb_ant': "merfish_202404091526/region_1",
    'P28DRb_ant': "merfish_202404091526/region_0",
    
    # exp 4-2
    'P21NRb_pos': "merfish_202404091601/region_0",
    'P21DRb_pos': "merfish_202404091601/region_2",
    'P28DRb_pos': "merfish_202404091601/region_1",    
}

merfish_datasets['P14NR_ant'] = merfish_datasets['P14NRa_ant']
merfish_datasets['P28NR_ant'] = merfish_datasets['P28NRa_ant']

merfish_datasets['P14NR_pos'] = merfish_datasets['P14NRa_pos']
merfish_datasets['P28NR_pos'] = merfish_datasets['P28NRa_pos']

merfish_datasets['P21NR_ant'] = merfish_datasets['P21NRa_ant']
merfish_datasets['P21DR_ant'] = merfish_datasets['P21DRa_ant']
merfish_datasets['P28DR_ant'] = merfish_datasets['P28DRa_ant']

merfish_datasets['P21NR_pos'] = merfish_datasets['P21NRa_pos']
merfish_datasets['P21DR_pos'] = merfish_datasets['P21DRa_pos']
merfish_datasets['P28DR_pos'] = merfish_datasets['P28DRa_pos']


# note that the index and the labels are pretty complicated but that is true
merfish_datasets_params = { 
    # 'P14NRa_ant': {'rotation':  90},
    # 'P14NRa_pos': {'rotation':  90},
    
#     'P14NRb_ant': {'rotation': -90},
#     'P14NRb_pos': {'rotation': -90},
    
#     'P21NRb_ant': {'rotation':  50},
#     'P21DRb_ant': {'rotation':  50},  
#     'P21NRb_pos': {'rotation':-130},
#     'P21DRb_pos': {'rotation':-130},  
    
    
    # right hemi + 45; left hemi -45
    'P28NRa_ant': {'rotation':  90, 'rotation2': -45},
    'P28NRa_pos': {'rotation': -90, 'rotation2':  45},
    
    'P28DRa_ant': {'rotation': -90, 'rotation2': -45},
    'P28DRa_pos': {'rotation': -90, 'rotation2': -45},    
    
    'P28NRb_ant': {'rotation':  90, 'rotation2': -45},
    'P28NRb_pos': {'rotation':  90, 'rotation2': -45},
    
    'P28DRb_ant': {'rotation': -30, 'rotation2':  60},
    'P28DRb_pos': {'rotation':-200, 'rotation2':  60},
}
