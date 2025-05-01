Setup:
  cdvars: ['co2_flx', 'camg_flx', 'totcat_flx', 'rockdiss']
  group_vars: ['dustrad', 'dustrate_ton_ha_yr']
  csv_sil: meanAnn_gbas_shortRun_base_lowFert_gs+apprate_v0.csv
  csv_cc: meanAnn_cc_shortRun_base_lowFert_gs+apprate_v0.csv
  multiyear_sil: False
  multiyear_cc: False

Removal accounting:
  time_horizon: 15.0
  cf_apprate_fixed: 0.6
  cf_dustrad_fixed: 100

Emissions (silicate):
  p80_input: 1300.0
  truck_km: 500.0
  barge_km: 0.0
  barge_diesel_km: 0
  Efactor_org: MRO

Emissions (calcite):
  p80_input: 1300.0
  truck_km: 500.0
  barge_km: 0.0
  barge_diesel_km: 0
  Efactor_org: MRO

