### dataset downloaded at the following link:
# https://zenodo.org/records/7142385
### downloaded 21/02/2024

### The dataset contains the rain gauge hourly rainfall series used in the
###  paper "How well does a convection-permitting climate model represent the
###  reverse orographic effect of extreme hourly precipitation?"

### Original dataset was stored in mat format. mat files were transformed in txt
###  format and then loaded in R. Each rain gauge series was saved in one R
###  variable, organized as a list L with five objects:
###   L$name: the identification name of the rain gauge station
###   L$vals_mm: series of hourly rainfall in millimeter
###   L$time_utc: time steps series, in UTC time
###   L$elev_m: elevation of the station, in m a.s.l.
###   L$xy_utm: station coordinates X and Y in meter in the Reference system WGS84/UTM zone 32N # nolint: line_length_linter.


# mat code mat2csv:

# listing = dir()
# listing = struct2cell(listing)
# listing = listing(1,3:176)

# for i = 1:174
#     name = listing{i}
#     load(name)
#     name = erase(name, ".mat")
#     vals_mm = S.vals_mm
#     fileName = strcat(name, "vals_mm")
#     writematrix(vals_mm, fileName)
#     time_utc = S.time_utc
#     fileName = strcat(name, "time_utc")
#     writematrix(time_utc, fileName)
#     xyUTM_elevM = [S.xy_utm, S.elev_m]
#     fileName = strcat(name, "xyUTM_elevM")
#     writematrix(xyUTM_elevM, fileName)
# end


# R code csv2R:

# #names = dir()
# #names = names[names!= "csv2RData.R" & names!= "mat2csv.m"]
# #save(names, file = "names.RData")
# load("names.RData")

# for (i in 1:length(names)) {
#   name = sub(".mat", "", names[i])
#   L = list()
#   vals_mm = read.csv(
#     file = paste(name, "vals_mm.txt", sep = ""), header = FALSE
#   )
#   vals_mm = vals_mm$V1
#   time_utc = read.csv(
#     file = paste(name, "time_utc.txt", sep = ""), header = FALSE
#   )
#   time_utc = time_utc$V1
#   xyUTM_elevM = read.csv(
#     file = paste(name, "xyUTM_elevM.txt", sep = ""), header = FALSE
#   )

#   L[["vals_mm"]] = vals_mm
#   L[["time_utc"]] = time_utc
#   L[["elev_m"]] = xyUTM_elevM$V3
#   L[["xy_utm"]] = c( xyUTM_elevM$V1, xyUTM_elevM$V2 )
#   L[["name"]] = name

#   save(L, file = paste(name, ".RData", sep = ""))
# }

# all(sub(".RData", "", dir("./rain_gauge_data_R")) == sub(".mat", "", names))

#' @import lubridate


### get daily data
# initialize object
daily_rainVTA_list = list()

# get data
for (i in dir("./rain_gauge_data_R")) {
  load(paste("./rain_gauge_data_R/", i, sep = ""))
  daily_rainVTA_list[[L$name]] = list()
  L$time_utc = lubridate::dmy_hms(L$time_utc)
  start_day = L$time_utc[1]
  end_day = tail(L$time_utc, 1)
  values = double(as.numeric(difftime(end_day, start_day), unit = "days") + 1)
  dates = character(length(values))
  origin_time = lubridate::as_datetime(lubridate::as_date(start_day))
  for (j in 1:length(values)) {
    values[j] = sum(L$vals_mm[
      L$time_utc >= origin_time + 3600 * 24 * (j - 1) &
        L$time_utc <= origin_time + 3600 * 24 * j - 1
    ])
    dates[j] = as.character(origin_time + 3600 * 24 * (j - 1))
  }
  dates = lubridate::as_date(dates)
  daily_rainVTA_list[[L$name]][["values_mm"]] = values
  daily_rainVTA_list[[L$name]][["dates_utc"]] = dates
  daily_rainVTA_list[[L$name]][["xyz"]] = c(L$xy_utm, L$elev_m)
}

### get overlapping time between stations
# initialize object
min_date = lubridate::ymd("0001-01-01")
max_date = lubridate::ymd("2999-12-31")

for (i in names(daily_rainVTA_list)) {
  if (daily_rainVTA_list[[i]]$dates_utc[1] > min_date) {
    min_date = daily_rainVTA_list[[i]]$dates_utc[1]
  }
  if (tail(daily_rainVTA_list[[i]]$dates_utc, 1) < max_date) {
    max_date = tail(daily_rainVTA_list[[i]]$dates_utc, 1)
  }
}
#min_date
#max_date

# get years
years = matrix(
  c(
    "2001-01-01", "2002-01-01", "2003-01-01", "2004-01-01",
    "2005-01-01", "2006-01-01", "2007-01-01", "2008-01-01",
    "2001-12-31", "2002-12-31", "2003-12-31", "2004-12-31",
    "2005-12-31", "2006-12-31", "2007-12-31", "2008-12-30"
  ), byrow = TRUE, nrow = 2
)

rainVTA_dataset = list()
rainVTA_dataset[["W"]] = list()
rainVTA_dataset[["S"]] = matrix(
  0,
  nrow = length(names(daily_rainVTA_list)),
  ncol = length(daily_rainVTA_list[[1]]$xyz)
)

for (i in 1:length(names(daily_rainVTA_list))) {
  rainVTA_dataset$W[[i]] = list()
  for (j in 1:ncol(years)) {
    rainVTA_dataset$W[[i]][[j]] = daily_rainVTA_list[[i]]$values_mm[
      daily_rainVTA_list[[i]]$dates_utc >= lubridate::as_date(years[1, j]) &
        daily_rainVTA_list[[i]]$dates_utc <= lubridate::as_date(years[2, j])
    ]
    rainVTA_dataset$W[[i]][[j]] = rainVTA_dataset$W[[i]][[j]][
      rainVTA_dataset$W[[i]][[j]] != 0
    ]
    if (length(rainVTA_dataset$W[[i]][[j]]) == 0) {
      rainVTA_dataset$W[[i]][[j]] = c(1e-3)
    }
  }
  rainVTA_dataset$S[i, ] = daily_rainVTA_list[[i]]$xyz
}
######## following (i, j) have only 0's (added a 0.001; before this,
#        minimum value was 0.1)
# "4 1"
# "54 1"
# "156 1"
# "160 1"
# "8 2"
# "41 2"
# "71 2"
# "25 3"
# "115 3"
# "27 4"
# "36 5"
# "43 6"
# "28 7"
# "37 7"
# "38 7"
# "12 8"
# "16 8"
# "17 8"
# "32 8"
# "40 8"

usethis::use_data(rainVTA_dataset)

# ### min function for POSIXct
# min_POSIXct = function(...) {
#   # convert POSIXct into seconds from origin
#   seconds = as.numeric(c(...))

#   # get minimum
#   min_seconds = min(seconds)

#   # convert minimum into POSIXct
#   min_POSIXct = lubridate::as_datetime(min_seconds)

#   # return the oldest date
#   return(min_POSIXct)
# }

# ### max function for POSIXct
# max_POSIXct = function(...) {
#   # convert POSIXct into seconds from origin
#   seconds = as.numeric(c(...))

#   # get maximum
#   max_seconds = max(seconds)

#   # convert maximum into POSIXct
#   max_POSIXct = lubridate::as_datetime(max_seconds)

#   # return the latest date
#   return(max_POSIXct)
# }

# provinces = vector("character", length = length(names))
# provinces[1:23] = "bolzano"
# provinces[24:44] = "trento"
# provinces[45] = "belluno" # VE_0003
# provinces[46] = "belluno" # VE_0009
# provinces[47] = "belluno" # VE_0011
# provinces[48] = "belluno" # VE_0015
# provinces[49] = "belluno" # VE_0019
# provinces[50] = "belluno" # VE_0021
# provinces[51] = "belluno" # VE_0022
# provinces[52] = "belluno" # VE_0025
# provinces[53] = "belluno" # VE_0037
# provinces[54] = "belluno" # VE_0047
# provinces[55] = "belluno" # VE_0048
# provinces[56] = "belluno" # VE_0050
# provinces[57] = "belluno" # VE_0053
# provinces[58] = "belluno" # VE_0055
# provinces[59] = "belluno" # VE_0056
# provinces[60] = "belluno" # VE_0058
# provinces[61] = "belluno" # VE_0059
# provinces[62] = "belluno" # VE_0061
# provinces[63] = "belluno" # VE_0067
# provinces[64] = "verona" # VE_0071
# provinces[65] = "vicenza" # VE_0072
# provinces[66] = "vicenza" # VE_0073
# provinces[67] = "vicenza" # VE_0076
# provinces[68] = "vicenza" # VE_0077
# provinces[69] = "vicenza" # VE_0079
# provinces[70] = "verona" # VE_0080
# provinces[71] = "vicenza" # VE_0081
# provinces[72] = "vicenza" # VE_0082
# provinces[73] = "vicenza" # VE_0083
# provinces[74] = "verona" # VE_0087
# provinces[75] = "vicenza" # VE_0088
# provinces[76] = "belluno" # VE_0091
# provinces[77] = "belluno" # VE_0092
# provinces[78] = "belluno" # VE_0093
# provinces[79] = "rovigo" # VE_0096
# provinces[80] = "rovigo" # VE_0098
# provinces[81] = "rovigo" # VE_0099
# provinces[82] = "treviso" # VE_0100
# provinces[83] = "rovigo" # VE_0101
# provinces[84] = "treviso" # VE_0102
# provinces[85] = "verona" # VE_0104
# provinces[86] = "vicenza" # VE_0105
# provinces[87] = "padova" # VE_0106
# provinces[88] = "verona" # VE_0108
# provinces[89] = "padova" # VE_0110
# provinces[90] = "padova" # VE_0111
# provinces[91] = "rovigo" # VE_0112
# provinces[92] = "rovigo" # VE_0113
# provinces[93] = "rovigo" # VE_0114
# provinces[94] = "rovigo" # VE_0115
# provinces[95] = "rovigo" # VE_0116
# provinces[96] = "verona" # VE_0117
# provinces[97] = "verona" # VE_0118
# provinces[98] = "verona" # VE_0119
# provinces[99] = "verona" # VE_0120
# provinces[100] = "rovigo" # VE_0121
# provinces[101] = "padova" # VE_0122
# provinces[102] = "verona" # VE_0123
# provinces[103] = "verona" # VE_0124
# provinces[104] = "verona" # VE_0125
# provinces[105] = "verona" # VE_0126
# provinces[106] = "verona" # VE_0127
# provinces[107] = "verona" # VE_0128
# provinces[108] = "verona" # VE_0129
# provinces[109] = "vicenza" # VE_0134
# provinces[110] = "vicenza" # VE_0135
# provinces[111] = "udine" # VE_0136  ##########################################
# provinces[112] = "vicenza" # VE_0137
# provinces[113] = "vicenza" # VE_0139
# provinces[114] = "vicenza" # VE_0140
# provinces[115] = "padova" # VE_0142
# provinces[116] = "vicenza" # VE_0144
# provinces[117] = "vicenza" # VE_0145
# provinces[118] = "vicenza" # VE_0147
# provinces[119] = "vicenza" # VE_0148
# provinces[120] = "vicenza" # VE_0149
# provinces[121] = "padova" # VE_0151
# provinces[122] = "padova" # VE_0152
# provinces[123] = "venezia" # VE_0159
# provinces[124] = "venezia" # VE_0160
# provinces[125] = "venezia" # VE_0163
# provinces[126] = "venezia" # VE_0165
# provinces[127] = "venezia" # VE_0166
# provinces[128] = "venezia" # VE_0167
# provinces[129] = "venezia" # VE_0168
# provinces[130] = "padova" # VE_0169
# provinces[131] = "padova" # VE_0170
# provinces[132] = "padova" # VE_0175
# provinces[133] = "padova" # VE_0177
# provinces[134] = "venezia" # VE_0178
# provinces[135] = "padova" # VE_0179
# provinces[136] = "padova" # VE_0182
# provinces[137] = "treviso" # VE_0183
# provinces[138] = "treviso" # VE_0184
# provinces[139] = "treviso" # VE_0185
# provinces[140] = "treviso" # VE_0186
# provinces[141] = "treviso" # VE_0187
# provinces[142] = "treviso" # VE_0188
# provinces[143] = "treviso" # VE_0189
# provinces[144] = "vicenza" # VE_0190
# provinces[145] = "vicenza" # VE_0191
# provinces[146] = "vicenza" # VE_0192
# provinces[147] = "treviso" # VE_0195
# provinces[148] = "treviso" # VE_0196
# provinces[149] = "treviso" # VE_0197
# provinces[150] = "belluno" # VE_0199
# provinces[151] = "belluno" # VE_0200
# provinces[152] = "trento" # VE_0203 ##########################################
# provinces[153] = "belluno" # VE_0216
# provinces[154] = "belluno" # VE_0217
# provinces[155] = "belluno" # VE_0219
# provinces[156] = "rovigo" # VE_0221
# provinces[157] = "belluno" # VE_0223
# provinces[158] = "belluno" # VE_0224
# provinces[159] = "treviso" # VE_0227
# provinces[160] = "venezia" # VE_0230
# provinces[161] = "rovigo" # VE_0231
# provinces[162] = "vicenza" # VE_0232
# provinces[163] = "padova" # VE_0234
# provinces[164] = "belluno" # VE_0235
# provinces[165] = "belluno" # VE_0236
# provinces[166] = "belluno" # VE_0237
# provinces[167] = "belluno" # VE_0238
# provinces[168] = "belluno" # VE_0239
# provinces[169] = "treviso" # VE_0240
# provinces[170] = "belluno" # VE_0246
# provinces[171] = "belluno" # VE_0247
# provinces[172] = "vicenza" # VE_0248
# provinces[173] = "verona" # VE_0251
# provinces[174] = "venezia" # VE_0252