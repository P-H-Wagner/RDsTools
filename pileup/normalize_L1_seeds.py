import ROOT
ROOT.EnableImplicitMT()

# 2018D_HLT_Mu7_IP4_run_from_320673_to_322349.root    358254225
# 2018D_HLT_Mu9_IP5_run_from_320673_to_322349.root    323236829
# 2018D_HLT_Mu9_IP6_run_from_320673_to_322349.root    130402521    
# 2018D_HLT_Mu12_IP6_run_from_320673_to_322349.root    62994349
# 
# 2018D_HLT_Mu7_IP4_run_from_322349_to_324295.root    196614520
# 2018D_HLT_Mu8_IP5_run_from_322349_to_324295.root    94550465
# 2018D_HLT_Mu9_IP5_run_from_322349_to_324295.root    70191581
# 2018D_HLT_Mu9_IP6_run_from_322349_to_324295.root    29475673
# 2018D_HLT_Mu12_IP6_run_from_322349_to_324295.root    27890368
# 
# 2018D_HLT_Mu7_IP4_run_from_324295_to_325175.root    239917306
# 2018D_HLT_Mu8_IP5_run_from_324295_to_325175.root    49071141
# 2018D_HLT_Mu9_IP5_run_from_324295_to_325175.root    32189579
# 2018D_HLT_Mu12_IP6_run_from_324295_to_325175.root    11732919

# 358254225+323236829+130402521+62994349+196614520+94550465+70191581+29475673+27890368+239917306+49071141+32189579+11732919

samples = dict()

samples['2018D'] = ROOT.RDataFrame('tree','/pnfs/psi.ch/cms/trivcat/store/user/manzoni/ParkingBPH2018D_part2_triggers_3Aug21_v1/trigger_tree_*.root')
samples['2018D'].Count().GetValue()
# 1626521476 events

samples['2018D'].Filter('320673<=run & run< 322349').Filter('HLT_Mu7_IP4_ps>0' )                                                                                    .Snapshot('tree', '2018D_HLT_Mu7_IP4_run_from_320673_to_322349.root' )
# 358254225 events

samples['2018D'].Filter('320673<=run & run< 322349').Filter('HLT_Mu9_IP5_ps>0' ).Filter('HLT_Mu7_IP4_ps==0')                                                        .Snapshot('tree', '2018D_HLT_Mu9_IP5_run_from_320673_to_322349.root' )
# 323236829 events

samples['2018D'].Filter('320673<=run & run< 322349').Filter('HLT_Mu9_IP6_ps>0' ).Filter('HLT_Mu7_IP4_ps==0').Filter('HLT_Mu9_IP5_ps==0')                            .Snapshot('tree', '2018D_HLT_Mu9_IP6_run_from_320673_to_322349.root' )
# 130402521 events

samples['2018D'].Filter('320673<=run & run< 322349').Filter('HLT_Mu12_IP6_ps>0').Filter('HLT_Mu7_IP4_ps==0').Filter('HLT_Mu9_IP5_ps==0').Filter('HLT_Mu9_IP6_ps==0').Snapshot('tree', '2018D_HLT_Mu12_IP6_run_from_320673_to_322349.root')
# 62994349 events

######################################

samples['2018D'].Filter('322349<=run & run< 324295').Filter('HLT_Mu7_IP4_ps>0' )                                                                                                                .Snapshot('tree', '2018D_HLT_Mu7_IP4_run_from_322349_to_324295.root' )
# 196614520 events

samples['2018D'].Filter('322349<=run & run< 324295').Filter('HLT_Mu8_IP5_ps>0' ).Filter('HLT_Mu7_IP4_ps==0')                                                                                    .Snapshot('tree', '2018D_HLT_Mu8_IP5_run_from_322349_to_324295.root' )
# 94550465 events

samples['2018D'].Filter('322349<=run & run< 324295').Filter('HLT_Mu9_IP5_ps>0' ).Filter('HLT_Mu7_IP4_ps==0').Filter('HLT_Mu8_IP5_ps==0')                                                        .Snapshot('tree', '2018D_HLT_Mu9_IP5_run_from_322349_to_324295.root' )
# 70191581 events

samples['2018D'].Filter('322349<=run & run< 324295').Filter('HLT_Mu9_IP6_ps>0' ).Filter('HLT_Mu7_IP4_ps==0').Filter('HLT_Mu8_IP5_ps==0').Filter('HLT_Mu9_IP5_ps==0')                            .Snapshot('tree', '2018D_HLT_Mu9_IP6_run_from_322349_to_324295.root' )
# 29475673 events

samples['2018D'].Filter('322349<=run & run< 324295').Filter('HLT_Mu12_IP6_ps>0').Filter('HLT_Mu7_IP4_ps==0').Filter('HLT_Mu8_IP5_ps==0').Filter('HLT_Mu9_IP5_ps==0').Filter('HLT_Mu9_IP6_ps==0').Snapshot('tree', '2018D_HLT_Mu12_IP6_run_from_322349_to_324295.root')
# 27890368 events

######################################

samples['2018D'].Filter('324295<=run & run<=325175').Filter('HLT_Mu7_IP4_ps>0' )                                                                                    .Snapshot('tree', '2018D_HLT_Mu7_IP4_run_from_324295_to_325175.root' )
# 239917306 events

samples['2018D'].Filter('324295<=run & run<=325175').Filter('HLT_Mu8_IP5_ps>0' ).Filter('HLT_Mu7_IP4_ps==0')                                                        .Snapshot('tree', '2018D_HLT_Mu8_IP5_run_from_324295_to_325175.root' )
# 49071141 events

samples['2018D'].Filter('324295<=run & run<=325175').Filter('HLT_Mu9_IP5_ps>0' ).Filter('HLT_Mu7_IP4_ps==0').Filter('HLT_Mu8_IP5_ps==0')                            .Snapshot('tree', '2018D_HLT_Mu9_IP5_run_from_324295_to_325175.root' )
# 32189579 events

samples['2018D'].Filter('324295<=run & run<=325175').Filter('HLT_Mu12_IP6_ps>0').Filter('HLT_Mu7_IP4_ps==0').Filter('HLT_Mu8_IP5_ps==0').Filter('HLT_Mu9_IP5_ps==0').Snapshot('tree', '2018D_HLT_Mu12_IP6_run_from_324295_to_325175.root')
# 11732919 events

