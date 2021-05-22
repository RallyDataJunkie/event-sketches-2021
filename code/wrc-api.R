library(jsonlite)
library(stringr)
library(dplyr)
library(purrr)

get_active_season = function(active_season_url="https://api.wrc.com/contel-page/83388/calendar/active-season/") {
  jsonlite::fromJSON(active_season_url)$rallyEvents$items
}

get_eventId_from_name = function(season, name){
  season[str_detect(season$name,
                    regex(name, ignore_case = T)), 'id']
}

results_api = 'https://api.wrc.com/results-api'

get_itinerary = function(eventId) {
  jsonlite::fromJSON(paste0(results_api,"/rally-event/",
                            eventId, "/itinerary"))$itineraryLegs
}

get_startlist = function(eventId, startListId) {
  startlist_url = paste0(results_api, '/rally-event/',
                         eventId,'/start-list-external/', startListId)
  
  startlist = jsonlite::fromJSON(startlist_url)$startListItems
  
  # Order the startlist dataframe by start order
  startlist %>% arrange(order)
}

get_startlist_id = function(itinerary, itinerarySectionId){
  sections = get_sections(itinerary)
  itineraryLegId = sections[sections$itinerarySectionId==itinerarySectionId,
                            'itineraryLegId']
  itinerary[itinerary$itineraryLegId==itineraryLegId,'startListId']
}

get_sections = function(itinerary){
  sections = do.call(rbind, itinerary$itinerarySections)
  sections
}

get_controls = function(sections){
  # Name each row in the list of dataframes we want to bind
  names(sections$controls) = sections$itinerarySectionId
  
  controls = sections$controls %>%
    # Ensure that we create an identifer column (uses list names)
    bind_rows(.id='itinerarySectionId')
  
  controls
}

get_stages = function(sections){
  # Name each row in the list of dataframes we want to bind
  names(sections$stages) = sections$itinerarySectionId
  
  stages = sections$stages %>%
    # Ensure that we create an identifer column (uses list names)
    bind_rows(.id='itinerarySectionId')
  
  stages
}

get_stage_list = function(stages){
  stage_list = stages$stageId
  stage_list
}

# https://stackoverflow.com/a/19265431/454773
get_stages_lookup = function(stages,
                             fromCol='code',  toCol='stageId'){
  stages_lookup = stages[[toCol]]
  names(stages_lookup) = stages[[fromCol]]
  stages_lookup
}

get_stage_id = function(stages, sname, typ='code'){
  # code, name
  if (typ=='code')
    stageId = stages[stages[typ] == sname, 'stageId']
  else
    stageId = stages[stringr::str_detect(stages[[typ]], sname), 'stageId']
  stageId
}

get_stage_info = function(stages, sid, typ='stageId', clean=TRUE){
  # stageId, code
  name=stages[stages[typ] == sid, 'name']
  distance=stages[stages[typ] == sid, 'distance']
  if (clean)
    stringr::str_replace(name, ' \\(Live TV\\)', '')
  
  c(name=name, distance=distance)
}

get_rally_entries = function(eventId) {
  cars_url = paste0(results_api, '/rally-event/',
                    eventId,'/cars')
  jsonlite::fromJSON(cars_url)
}

get_drivers = function(entries){
  drivers = do.call(cbind, entries$driver)
  drivers
}

get_codrivers = function(entries){
  codrivers = bind_cols(entries$codriver)
  codrivers
}

get_person_id = function(persons, sname,
                         typ='fullName', ret='personId'){
  
  if (!(typ %in% colnames(persons)))
    typ='driverfullname'
  
  if (!(ret %in% colnames(persons)))
    ret='driverId'
  
  # code, fullName
  if (typ=='code')
    personId = persons[persons[typ]==sname, ret]
  else
    personId = persons[str_detect(persons[[typ]],
                                  regex(sname,
                                        ignore_case = T)),
                       ret]
  personId
}

get_car_data = function(entries){
  cols = c('entryId', 'driverId', 'codriverId', 'manufacturerId',
           'vehicleModel','eligibility', 'classname','manufacturer',
           'entrantname', 'groupname', 
           'identifier', 'drivername', 'code',
           'driverfullname', 'codrivername','codriverfullname'
  )
  entries = entries %>%
    rowwise() %>% 
    mutate(classname = eventClasses$name) %>%
    mutate(manufacturer = manufacturer$name) %>%
    mutate(entrantname = entrant$name) %>%
    mutate(groupname = group$name) %>%
    mutate(drivername = driver$abbvName) %>%
    mutate(driverfullname = driver$fullName) %>%
    mutate(codrivername = codriver$abbvName) %>%
    mutate(codriverfullname = codriver$fullName) %>%
    mutate(code = driver$code) %>%
    select(all_of(cols))
  
  # If we don't cast, it's a non-rankable rowwise df
  as.data.frame(entries)
}

get_penalties = function(eventId) {
  penalties_url = paste0(results_api, '/rally-event/',
                         eventId, '/penalties')
  jsonlite::fromJSON(penalties_url)
}

get_retirements = function(eventId) {
  retirements_url = paste0(results_api, '/rally-event/',
                           eventId, '/retirements')
  jsonlite::fromJSON(retirements_url)
}


get_result = function(eventId) {
  result_url = paste0(results_api, '/rally-event/',
                      eventId,'/result')
  
  jsonlite::fromJSON(result_url)
}


get_stage_winners = function(eventId) {
  stage_winners_url = paste0(results_api, '/rally-event/',
                             eventId,'/stage-winners')
  
  jsonlite::fromJSON(stage_winners_url)
}

get_overall_result = function(eventId, stageId) {
  overall_url = paste0(results_api, '/rally-event/',
                       eventId, '/stage-result/stage-external/',
                       stageId)
  jsonlite::fromJSON(overall_url) %>%
    # Also add in the stage ID
    mutate(stageId = stageId)
}

get_overall_result2 = function(stageId, eventId) {
  get_overall_result(eventId, stageId)
}

get_multi_overall = function(stage_list){
  multi_overall = stage_list %>%
    map(get_overall_result2, eventId=eventId) %>% 
    bind_rows()
  multi_overall
}

get_stage_times = function(eventId, stageId) {
  stage_times_url = paste0(results_api, '/rally-event/',
                           eventId, '/stage-times/stage-external/',
                           stageId)
  jsonlite::fromJSON(stage_times_url)
}

get_stage_times2 = function(stageId, eventId) {
  get_stage_times(eventId, stageId)
}

get_multi_stage_times = function(stage_list){
  multi_stage_times = stage_list %>%
    map(get_stage_times2, eventId=eventId) %>% 
    bind_rows()
  multi_stage_times
}

get_multi_stage_times_wide = function(multi_stage_times, stage_list){
  stage_times_cols = c('entryId', 'stageId', 'elapsedDurationMs')
  
  multi_stage_times_wide = multi_stage_times %>% 
    select(all_of(stage_times_cols)) %>%
    mutate(elapsedDurationS = elapsedDurationMs / 1000) %>%
    select(-elapsedDurationMs) %>%
    group_by(entryId) %>%
    tidyr::spread(key = stageId,
                  value = elapsedDurationS) %>%
    select(c('entryId', as.character(stage_list))) %>%
    # If we don't cast, it's a
    # non-rankable rowwise df
    as.data.frame()
  
  multi_stage_times_wide
}

get_multi_stage_positions_wide = function(multi_stage_times, stage_list){
  stage_positions_cols = c('entryId', 'stageId', 'position')
  
  multi_stage_positions_wide = multi_stage_times %>% 
    select(all_of(stage_positions_cols)) %>%
    group_by(entryId) %>%
    tidyr::spread(key = stageId,
                  value = position) %>%
    select(c('entryId', as.character(stage_list))) %>%
    # If we don't cast, it's a
    # non-rankable rowwise df
    as.data.frame()
}

get_multi_stage_generic_wide = function(multi_stage_generic, stage_list,
                                        wide_val, group_key='entryId',
                                        spread_key='stageId'){
  
  stage_times_cols = c(group_key, spread_key, wide_val )
  
  if (wide_val=='elapsedDurationMs') {
    multi_stage_times_wide = multi_stage_times %>% 
      select(all_of(stage_times_cols)) %>%
      mutate(elapsedDurationS = elapsedDurationMs / 1000) %>%
      select(-elapsedDurationMs)
    wide_val = 'elapsedDurationS'
  }
  
  multi_stage_generic_wide = multi_stage_generic %>% 
    select(all_of(stage_times_cols)) %>%
    # group_by_at lets us pass in the grouping column by variable
    group_by_at(group_key) %>%
    tidyr::spread(key = spread_key,
                  value = wide_val) %>%
    select( c(group_key, as.character(stage_list))) %>%
    # If we don't cast, it's a
    # non-rankable rowwise df
    as.data.frame()
  
  multi_stage_generic_wide
}


get_splits = function(eventId, stageId){
  splits_url=paste0(results_api, '/rally-event/', eventId,
                    '/split-times/stage-external/', stageId)
  
  jsonlite::fromJSON(splits_url)
}

get_driver_splits = function(splits){
  driver_splits = splits$entrySplitPointTimes$splitPointTimes %>%
    bind_rows() %>%
    mutate(elapsedDurationS = elapsedDurationMs / 1000) %>%
    select(-elapsedDurationMs)
  driver_splits
}

get_split_cols = function(splits){
  split_cols =  as.character(arrange(splits$splitPoints, distance)$splitPointId)
  split_cols
}

get_driver_splits_wide = function(driver_splits, splits){
  split_cols =  get_split_cols(splits)
  splits_cols = c('entryId', 'splitPointId', 'elapsedDurationS')
  
  driver_splits_wide = driver_splits %>% 
    group_by(entryId) %>%
    select(all_of(splits_cols)) %>%
    tidyr::spread(key = splitPointId,
                  value = elapsedDurationS) %>%
    select(all_of(c('entryId', split_cols))) %>%
    # If we don't cast, it's a
    # non-rankable rowwise df
    as.data.frame()
  driver_splits_wide
}


rebase = function(df, id, rebase_cols,
                  id_col='entryId', base=FALSE,  base_id=FALSE) {
  
  df_ =  df
  
  rebase_cols = as.character(rebase_cols)
  
  # The rebase values are the ones
  # we want to subtract from each row
  rebase_vals = c(df[df[[id_col]]==id, rebase_cols])
  
  # Do the rebasing
  df_[,rebase_cols] =  df[,rebase_cols] - rebase_vals
  
  df_[[id_col]] = df[[id_col]]
  
  # Return just the rebased and identifier columns or the
  # whole dataframe
  cols = rebase_cols
  if (base_id)
    cols = c(id_col, cols)
  if (base)
    df_ %>% select(cols)
  else
    df_
}