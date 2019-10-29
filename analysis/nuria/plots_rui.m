
%-- 10/24/2019 1:06 PM --%
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI12\D5\BMI_online190930T152419.mat')
sum(data.selfHits)
Holo15total = [18];
Holo15hpm = [18/40];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI13\D5\BMI_online190930T180325.mat')
sum(data.selfHits)
Holo15total = [Holo15total, 69];
Holo15hpm = [Holo15hpm, 69/40];
plot(smoothdata(data.selfHits,'movmean',10000))
plot(data.selfHits)
plot(smoothdata(data.selfHits,'movmean',10000))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI12\D5\BMI_online190930T152419.mat')
plot(smoothdata(data.selfHits,'movmean',10000))
kk = data.selfHits;
selfhits = data.selfHits;
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI13\D5\BMI_online190930T180325.mat')
selfhits = [selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI16\D5\BMI_online190930T202917.mat')
sum(data.selfHits)
Holo15total = [Holo15total, 19];
Holo15hpm = [Holo15hpm, 19/40];
selfhits = [selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI12\D6\BMI_online191001T152524.mat')
sum(data.selfHits)
Holo15total = [Holo15total, 35];
Holo15hpm = [Holo15hpm, 35/40];
selfhits = [selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI13\D6\BMI_online191001T180737.mat')
sum(data.selfHits)
Holo15total = [Holo15total, 20];
Holo15hpm = [Holo15hpm, 20/40];
selfhits = [selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_online191001T203900.mat')
sum(data.selfHits)
plot(smoothdata(data.selfHits,'movmean',10000))
Holo15total = [Holo15total, 42];
selfhits = [selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_online191001T212200.mat')
sum(data.selfHits)
Holo15hpm = [Holo15hpm, 85/80];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191003\NVI12\D8\BMI_online191003T153014.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191003\NVI13\D8\BMI_online191003T175350.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191003\NVI13\D8\BMI_online191003T184251.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191003\NVI13\D8\BMI_online191003T175350.mat')
Holo15total = [Holo15total, 19];
selfhits = [selfhits; data.selfHits];
Holo15hpm = [Holo15hpm, 50/80];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI12\D10\BMI_online191005T155340.mat')
sum(data.selfHits)
Holo15total = [Holo15total, 29];
selfhits = [selfhits; data.selfHits];
Holo15hpm = [Holo15hpm, 29/40];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI13\D10\BMI_online191005T183253.mat')
sum(data.selfHits)
Holo15total = [Holo15total, 25];
selfhits = [selfhits; data.selfHits];
Holo15hpm = [Holo15hpm, 25/40];
clear all
clc
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI12\D5\BMI_online190930T152419.mat')
sum(data.selfHits)
nvi12_holo_15_total = [18];
nvi12_holo_15_jpm = [18/40];
selfhits = [data.selfHits];
nvi12_selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI13\D5\BMI_online190930T180325.mat')
nvi13_holo_15_total = [69];
sum(data.selfHits)
nvi13_holo_15_hpm = [69/40];
nvi13_selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI16\D5\BMI_online190930T202917.mat')
sum(data.selfHits)
nvi16_holo_15_total = [19];
nvi16_holo_15_hpm = [19/40];
nvi16_selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI12\D6\BMI_online191001T152524.mat')
sum(data.selfHits)
nvi12_holo_15_total = [nvi12_holo_15_total, 35];
nvi12_holo_15_hpm = [nvi12_holo_15_hpm, 35/40];
nvi12_selfhits = [nvi12_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI13\D6\BMI_online191001T180737.mat')
sum(data.selfHits)
nvi12_holo_15_total = [nvi12_holo_15_total, 35];
nvi12_holo_15_hpm = [nvi12_holo_15_hpm, 35/40];
nvi12_selfhits = [nvi12_selfhits; data.selfHits];
nvi12_holo_15_hpm[end] = []
nvi12_holo_15_hpm(end) = []
nvi12_holo_15_total(end) = []
nvi12_selfhits(end,:) = [];
sum(data.selfHits)
nvi13_holo_15_total = [nvi13_holo_15_total, 20];
nvi13_holo_15_hpm = [nvi13_holo_15_hpm, 20/40];
nvi13_selfhits = [nvi13_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_online191001T203900.mat')
sum(data.selfHits)
nvi16_holo_15_total = [nvi16_holo_15_total, 42];
nvi16_holo_15_hpm = [nvi16_holo_15_hpm, 85/80];
nvi16_selfhits = [nvi16_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191003\NVI13\D8\BMI_online191003T175350.mat')
sum(data.selfHits)
nvi13_holo_15_total = [nvi13_holo_15_total, 19];
nvi13_selfhits = [nvi13_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191003\NVI13\D8\BMI_online191003T184251.mat')
sum(data.selfHits)
nvi13_holo_15_hpm = [nvi13_holo_15_hpm, 50/80];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI12\D10\BMI_online191005T155340.mat')
sum(data.selfHits)
nvi12_holo_15_total = [nvi12_holo_15_total, 29];
nvi12_holo_15_hpm = [nvi12_holo_15_hpm, 29/40];
nvi12_selfhits = [nvi12_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI13\D10\BMI_online191005T183253.mat')
sum(data.selfHits)
nvi13_holo_15_total = [nvi13_holo_15_total, 25];
nvi13_holo_15_hpm = [nvi13_holo_15_hpm, 25/40];
nvi13_selfhits = [nvi13_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI16\D10\BMI_online191005T213208.mat')
sum(data.selfHits)
nvi16_holo_15_total = [nvi16_holo_15_total, 37];
nvi16_selfhits = [nvi16_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI16\D10\BMI_online191005T221448.mat')
sum(data.selfHits)
nvi16_holo_15_hpm = [nvi16_holo_15_hpm, 77/80];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI12\D12\BMI_online191007T153922.mat')
sum(data.selfHits)
nvi12_holo_15_total = [nvi12_holo_15_total, 22];
nvi12_holo_15_hpm = [nvi12_holo_15_hpm, 22/40];
nvi12_selfhits = [nvi12_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI12\D12\BMI_online191007T155538.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI13\D12\BMI_online191007T202814.mat')
sum(data.selfHits)
nvi13_holo_15_total = [nvi13_holo_15_total, 6];
nvi13_holo_15_hpm = [nvi13_holo_15_hpm, 6/40];
nvi13_selfhits = [nvi13_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI16\D12\BMI_online191007T233928.mat')
sum(data.selfHits)
nvi16_holo_15_total = [nvi16_holo_15_total, 15];
nvi16_holo_15_hpm = [nvi16_holo_15_hpm, 15/40];
nvi16_selfhits = [nvi16_selfhits; data.selfHits];
[2 99]
y
n
[25 27 28]
n
y
n
load('H:\holobmi_H\191024\NVI12\D24\holostim_seq191024T151618.mat')
load('H:\holobmi_H\191024\NVI12\D24\holostim_seq191024T151920.mat')
uiopen('H:\holobmi_H\191024\NVI12\D24\im\holostim_seq\holostim_seq_191024T151755-070\holostim_seq_191024T151755-070_Cycle00001_VoltageRecording_001.csv',1)
[task_settings] = define_BMI_task_settings(fb_bool);
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI12\D14\BMI_online191009T153223.mat')
sum(data.selfHits)
nvi12_holo_15_total = [nvi12_holo_15_total, 75];
nvi12_holo_15_hpm = [nvi12_holo_15_hpm, 75/40];
nvi12_selfhits = [nvi12_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI13\D14\BMI_online191009T192110.mat')
sum(data.selfHits)
nvi13_holo_15_total = [nvi13_holo_15_total, 22];
nvi13_holo_15_hpm = [nvi13_holo_15_hpm, 22/40];
nvi13_selfhits = [nvi13_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI16\D14\BMI_online191009T220103.mat')
sum(data.selfHits)
nvi16_holo_15_total = [nvi16_holo_15_total, 13];
nvi16_holo_15_hpm = [nvi16_holo_15_hpm, 13/40];
nvi16_selfhits = [nvi16_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI12\D16\BMI_online191011T151240.mat')
sum(data.selfHits)
nvi12_holo_15_total = [nvi12_holo_15_total, 39];
nvi12_holo_15_hpm = [nvi12_holo_15_hpm, 39/40];
nvi12_selfhits = [nvi12_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI13\D16\BMI_online191011T175050.mat')
sum(data.selfHits)
nvi13_holo_15_total = [nvi13_holo_15_total, 10];
nvi13_holo_15_hpm = [nvi13_holo_15_hpm, 10/40];
nvi13_selfhits = [nvi13_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI16\D16\BMI_online191011T205824.mat')
sum(data.selfHits)
nvi16_holo_15_total = [nvi16_holo_15_total, 19];
nvi16_selfhits = [nvi16_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI16\D16\BMI_online191011T214141.mat')
sum(data.selfHits)
nvi16_holo_15_hpm = [nvi16_holo_15_hpm, 46/80];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191013\NVI12\D18\BMI_online191013T150710.mat')
sum(data.selfHits)
nvi12_holo_15_total = [nvi12_holo_15_total, 22];
nvi12_holo_15_hpm = [nvi12_holo_15_hpm, 22/40];
nvi12_selfhits = [nvi12_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191013\NVI13\D18\BMI_online191013T180942.mat')
sum(data.selfHits)
nvi13_holo_15_total = [nvi13_holo_15_total, 12];
nvi13_holo_15_hpm = [nvi13_holo_15_hpm, 12/40];
nvi13_selfhits = [nvi13_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191013\NVI16\D18\BMI_online191013T212923.mat')
sum(data.selfHits)
nvi16_holo_15_total = [nvi16_holo_15_total, 13];
nvi16_selfhits = [nvi16_selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191013\NVI16\D18\BMI_online191013T222954.mat')
sum(data.selfHits)
nvi16_holo_15_hpm = [nvi16_holo_15_hpm, 13/40];
nvi12_holo_15_total(end) = []
nvi12_holo_15_hpm(end) = []
nvi12_selfhits(end,:) = [];
nvi16_holo_15_total(end) = []
nvi16_holo_15_hpm(end) = []
nvi16_selfhits(end,:) = [];
run('Z:\Vivek\VivekNuria\HoloBMI_092319\191015\NVI12\D20\mainProt_single_channel.m')
nanmean(nvi12_holo_15_total)
nanmean(nvi13_holo_15_total)
nanmean(nvi16_holo_15_total)
vectorVTA = vectorHolo
vectorHolo = []
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191015\NVI12\D20\BMI_online191015T143410.mat')
sum(data.selfHits)
nvi12_holo_40_total = [12];
nvi12_holo_40_hpm = [12/40];
nvi12_40selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191015\NVI13\D20\BMI_online191015T175005.mat')
sum(data.selfHits)
nvi13_holo_40_total = [30];
nvi13_holo_40_hpm = [30/40];
nvi13_40selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191015\NVI16\D20\BMI_online191015T204909.mat')
sum(data.selfHits)
nvi16_holo_40_total = [7];
nvi16_holo_40_hpm = [7/40];
nvi16_40selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191017\NVI13\D22\BMI_online191017T184758.mat')
sum(data.selfHits)
nvi13_holo_40_total = [nvi13_holo_40_total, 8];
nvi13_holo_40_hpm = [nvi13_holo_40_hpm, 8/40];
nvi13_40selfhits = [nvi13_40selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191017\NVI16\D22\BMI_online191017T221038.mat')
sum(data.selfHits)
nvi16_holo_40_total = [nvi16_holo_40_total, 14];
nvi16_holo_40_hpm = [nvi16_holo_40_hpm, 14/40];
nvi16_40selfhits = [nvi16_40selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191006\NVI12\D11\BMI_online191006T152205.mat')
sum(data.selfHits)
nvi12_rr_15_total = [7]
nvi12_rr_15_hpm = [7/40]
nvi12_rrselfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191006\NVI13\D11\BMI_online191006T185125.mat')
sum(data.selfHits)
nvi13_rr_15_total = [12]
nvi13_rr_15_hpm = [12/40]
nvi13_rrselfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191006\NVI16\D11\BMI_online191006T212453.mat')
sum(data.selfHits)
nvi16_rr_15_total = [11]
nvi16_rr_15_hpm = [11/40]
nvi16_rrselfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191010\NVI12\D15\BMI_online191010T145643.mat')
sum(data.selfHits)
nvi12_rr_15_total = [nvi12_rr_15_total, 25]
nvi12_rr_15_hpm = [nvi12_rr_15_hpm, 25/40]
nvi12_rrselfhits = [nvi12_rrselfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191010\NVI13\D15\BMI_online191010T173433.mat')
sum(data.selfHits)
nvi13_rr_15_total = [nvi13_rr_15_total, 34]
nvi13_rr_15_hpm = [nvi13_rr_15_hpm, 34/40]
nvi13_rrselfhits = [nvi13_rrselfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191010\NVI16\D15\BMI_online191010T201719.mat')
sum(data.selfHits)
nvi16_rr_15_total = [nvi16_rr_15_total, 31]
nvi16_rr_15_hpm = [nvi16_rr_15_hpm, 31/40]
nvi16_rrselfhits = [nvi16_rrselfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191014\NVI12\D19\BMI_online191014T160227.mat')
sum(data.selfHits)
nvi12_rr_15_total = [nvi12_rr_15_total, 14]
nvi12_rr_15_hpm = [nvi12_rr_15_hpm, 14/40]
nvi12_rrselfhits = [nvi12_rrselfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191014\NVI13\D19\BMI_online191014T184703.mat')
sum(data.selfHits)
nvi13_rr_15_total = [nvi13_rr_15_total, 12]
nvi13_rr_15_hpm = [nvi13_rr_15_hpm, 12/40]
nvi13_rrselfhits = [nvi13_rrselfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191014\NVI16\D19\BMI_online191014T225326.mat')
sum(data.selfHits)
nvi16_rr_15_total = [nvi16_rr_15_total, 16]
nvi16_rr_15_hpm = [nvi16_rr_15_hpm, 16/40]
nvi16_rrselfhits = [nvi16_rrselfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191018\NVI12\D23\BMI_online191018T152337.mat')
sum(data.selfHits)
nvi12_rr_40_total = [31]
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191018\NVI12\D23\BMI_online191018T160622.mat')
sum(data.selfHits)
mean(diff(data.selfHits))
diff(data.selfHits)
clc
holo15total = [nanmean(nvi12_holo_15_total), nanmean(nvi13_holo_15_total), nanmean(nvi16_holo_15_total)]
holo15hpm = [nanmean(nvi12_holo_15_hpm), nanmean(nvi13_holo_15_hpm), nanmean(nvi16_holo_15_hpm)]
holo40total = [nanmean(nvi12_holo_40_total), nanmean(nvi13_holo_40_total), nanmean(nvi16_holo_40_total)]
holo15total = [nanmean(nvi12_rr_15_total), nanmean(nvi13_rr_15_total), nanmean(nvi16_rr_15_total)]
holo15hpm = [nanmean(nvi12_holo_15_hpm), nanmean(nvi13_holo_15_hpm), nanmean(nvi16_holo_15_hpm)]
holo15total = [nanmean(nvi12_holo_15_total), nanmean(nvi13_holo_15_total), nanmean(nvi16_holo_15_total)]
rr15total = [nanmean(nvi12_rr_15_total), nanmean(nvi13_rr_15_total), nanmean(nvi16_rr_15_total)]
rr15hpm = [nanmean(nvi12_rr_15_hpm), nanmean(nvi13_rr_15_hpm), nanmean(nvi16_rr_15_hpm)]
yy = [nanmean(nvi12_holo_15_total), nanmean(nvi12_rr_15_total), nanmean(nvi13_holo_15_total), nanmean(nvi13_rr_15_total), nanmean(nvi16_holo_15_total), nanmean(nvi16_rr_15_total)
];
xx = [1,2, 4,5, 7,8];
yy = [std(nvi12_holo_15_total)/sqrt(length(nvi12_holo_15_total)), std(nvi12_rr_15_total)/sqrt(length(nvi12_rr_15_total))), std(nvi13_holo_15_total)/sqrt(length(nvi13_holo_15_total)), std(nvi13_rr_15_total)/sqrt(length(nvi13_rr_15_total)), std(nvi16_holo_15_total)/sqrt(length(nvi16_holo_15_total)), std(nvi16_rr_15_total)/sqrt(length(nvi16_rr_15_total))]
yy = [std(nvi12_holo_15_total)/sqrt(length(nvi12_holo_15_total)), std(nvi12_rr_15_total)/sqrt(length(nvi12_rr_15_total)), std(nvi13_holo_15_total)/sqrt(length(nvi13_holo_15_total)), std(nvi13_rr_15_total)/sqrt(length(nvi13_rr_15_total)), std(nvi16_holo_15_total)/sqrt(length(nvi16_holo_15_total)), std(nvi16_rr_15_total)/sqrt(length(nvi16_rr_15_total))]
yy = [nanmean(nvi12_holo_15_total), nanmean(nvi12_rr_15_total), nanmean(nvi13_holo_15_total), nanmean(nvi13_rr_15_total), nanmean(nvi16_holo_15_total), nanmean(nvi16_rr_15_total)
]
err = [std(nvi12_holo_15_total)/sqrt(length(nvi12_holo_15_total)), std(nvi12_rr_15_total)/sqrt(length(nvi12_rr_15_total)), std(nvi13_holo_15_total)/sqrt(length(nvi13_holo_15_total)), std(nvi13_rr_15_total)/sqrt(length(nvi13_rr_15_total)), std(nvi16_holo_15_total)/sqrt(length(nvi16_holo_15_total)), std(nvi16_rr_15_total)/sqrt(length(nvi16_rr_15_total))]
errorbar(x,y,err)
errorbar(xx,yy,err)
bar(xx,yy)
hold on
errorbar(xx,yy,err)
xticklabels({'12H', '12rr', '13H', '13rr', '16H', '16rr'})
ttest(nanmean(nvi12_holo_15_total, nanmean(nvi12_rr_15_total)
ttest(nanmean(nvi12_holo_15_total, nanmean(nvi12_rr_15_total)))
ttest(nvi12_holo_15_total, nvi12_rr_15_total)
ttest2(nvi12_holo_15_total, nvi12_rr_15_total)
[h,p] = ttest2(nvi12_holo_15_total, nvi12_rr_15_total)
nvi12_holo_15_total
nvi12_rr_15_total
[h,p] = ttest2(nvi12_holo_15_hpm, nvi12_rr_15_hpm)
figure
plot(nvi12_holo_15_hpm)
plot(nvi12_holo_15_hpm + nvi16_holo_15_hpm)
figure
plot(nvi16_holo_15_hpm)
plot(nvi16_holo_15_total)
vline(20)
hline(20)
bar(xx,yy)
errorbar(xx,yy,err)
xticklabels({'12H', '12rr', '13H', '13rr', '16H', '16rr'})
hline(20)
bar(xx,yy)
errorbar(xx,yy,err)
hline(20)
bar(xx,yy)
hold on
errorbar(xx,yy,err)
hline(20)
xticklabels({'12H', '12rr', '13H', '13rr', '16H', '16rr'})
xticklabels({'12_H', '12_rr', '13_H', '13_rr', '16_H', '16_rr'})
xticklabels({'12_H', '12_C', '13_H', '13_C', '16_H', '16_C'})
xticklabels({'1_H', '1_C', '2_H', '2_C', '3_H', '3_C'})
plot(nvi12_rr_15_hpm + nvi16_rr_15_hpm)
plot(nvi12_rr_15_total + nvi16_rr_15_total)
plot(nvi12_rr_15_total + nvi13_rr_15_total + nvi16_rr_15_total)
load('H:\holobmi_H\191024\NVI12\D24\BMI_online191024T171438.mat')
sum(data.selfHits)
plot(smoothdata(data.selfHits,'movmean',10000))
[2 99]
y
n
[]
n
y
n
load('E:\holobmi_E\191024\NVI13\D24\holostim_seq191024T184053.mat')
uiopen('E:\holobmi_E\191024\NVI13\D24\im\holostim_seq\holostim_seq_191024T183938-076\holostim_seq_191024T183938-076_Cycle00001_VoltageRecording_001.csv',1)
vectorVTA = vectorHolo
vectorHolo = []
clear all
load('I:\Vivek\Imaging_data\191024\NVI12\D24\BMI_online191024T171438.mat')
sum(data.selfHits)
load('I:\Vivek\Imaging_data\191024\NVI12\D24\BMI_online191024T163039.mat')
sum(data.selfHits)
plot(smoothdata(data.selfHits,'movmean',10000))
kk = data.selfHits;
kk(find(kk==1))
find(kk==1)
diff(find(kk==1))
for i=1:length(kk)
end
hits = find(kk==1);
hits((diff(hits)==1)+1) = []
hits
hits = find(kk==1);
diff(hits)
hits((diff(hits)==1))
kk((diff(hits)==1))=0;
sum(kk)
sum(data.selfHits)
kk((diff(hits)==1))
kk = data.selfHits;
kk((diff(hits)==1))
hits = find(kk==1);
hits
kk((diff(hits)==1))
diff(hits)
hits(diff(hits)==1)
kk(hits(diff(hits)==1))=0
sum(kk)
plot(kk)
kk(hits(diff(hits)==1))=0;
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('E:\holobmi_E\191024\NVI13\D24\BMI_online191024T203705.mat')
sum(data.selfHits)
plot(smoothdata(data.selfHits,'movmean',10000))
clear all
close all
load('E:\holobmi_E\191024\NVI13\D24\BMI_online191024T203705.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('E:\holobmi_E\191024\NVI13\D24\BMI_online191024T195246.mat')
kk = data.selfHits;
clear all
%-- 10/24/2019 8:46 PM --%
load('E:\holobmi_E\191024\NVI13\D24\BMI_online191024T195246.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)==1))=0;
sum(kk)
load('E:\holobmi_E\191024\NVI13\D24\BMI_online191024T203705.mat')
sum(data.selfHits)
[2 99]
y
n
[]
n
y
n
load('H:\holobmi_H\191024\NVI16\D24\holostim_seq191024T211552.mat')
uiopen('H:\holobmi_H\191024\NVI16\D24\im\holostim_seq\holostim_seq_191024T211430-082\holostim_seq_191024T211430-082_Cycle00001_VoltageRecording_001.csv',1)
vectorVTA = vectorHolo
vectorHolo = []
load('H:\holobmi_H\191024\NVI16\D24\BMI_online191024T230420.mat')
sum(data.selfHits)
load('H:\holobmi_H\191024\NVI16\D24\BMI_online191024T222011.mat')
kk = data.selfHits;
%-- 10/24/2019 11:11 PM --%
load('H:\holobmi_H\191024\NVI16\D24\BMI_online191024T222011.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)==1))=0;
sum(kk)
plot(kk)
for i=1:length(kk)
if kk(i)==1
kk(i:5*30)=0
end
end
kk(i:5*30)
i=1
kk(i:5*30)
kk(i:5*30)=1
kk(i:5*30)
kk(i+1:5*30)=0
kk(i+1:5*30)
kk(i)
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)==1))=0;
for i=1:length(kk)
if kk(i)==1
kk(i+1:5*30)=0;
end
sum(kk)
kk = data.selfHits;
kk(i+1:5*30)=0;
kk = data.selfHits;
kk(hits(diff(hits)==1))=0;
hits = find(kk==1);
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=5*30))=0;
sum(kk)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191024\NVI12\D24\BMI_online191024T163039.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=5*30))=0;
sum(kk)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191024\NVI13\D24\BMI_online191024T195246.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=5*30))=0;
sum(kk)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191024\NVI12\D24\BMI_online191024T163039.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=2*30))=0;
sum(kk)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191024\NVI12\D24\BMI_online191024T171438.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=2*30))=0;
sum(kk)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191024\NVI13\D24\BMI_online191024T195246.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=2*30))=0;
sum(kk)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191024\NVI13\D24\BMI_online191024T203705.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=2*30))=0;
sum(kk)
load('H:\holobmi_H\191024\NVI16\D24\BMI_online191024T222011.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=2*30))=0;
sum(kk)
load('H:\holobmi_H\191024\NVI16\D24\BMI_online191024T230420.mat')
kk = data.selfHits;
hits = find(kk==1);
kk(hits(diff(hits)<=2*30))=0;
sum(kk)
clear all
close all
[2 99]
y
n
[]
n
y
n
load('E:\holobmi_E\191025\NVI12\D25\holostim_seq191025T132909.mat')
uiopen('E:\holobmi_E\191025\NVI12\D25\im\holostim_seq\holostim_seq_191025T132800-087\holostim_seq_191025T132800-087_Cycle00001_VoltageRecording_001.csv',1)
%-- 10/25/2019 2:43 PM --%
load('E:\holobmi_E\191025\NVI12\D25\BMI_online191025T143753.mat')
sum(data.holoHits)
load('E:\holobmi_E\191025\NVI12\D25\BMI_online191025T154101.mat')
sum(data.selfHits)
clear all
load('G:\VivekNuria\Code\HoloBMI\analysis\nuria\checking_itworks.mat')
[2 99]
y
n
[1]
n
y
n
load('H:\holobmi_H\191025\NVI13\D25\holostim_seq191025T162423.mat')
uiopen('H:\holobmi_H\191025\NVI13\D25\im\holostim_seq\holostim_seq_191025T162324-094\holostim_seq_191025T162324-094_Cycle00001_VoltageRecording_001.csv',1)
%-- 10/25/2019 4:40 PM --%
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191012\NVI12\D17\BMI_online191012T144049.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191012\NVI13\D17\BMI_online191012T170122.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191012\NVI16\D17\BMI_online191012T191132.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191008\NVI12\D13\BMI_online191008T152023.mat')
sum(data.selfHits)
19+39
58/2
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191008\NVI13\D13\BMI_online191008T175827.mat')
sum(data.selfHits)
48/2
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191008\NVI16\D13\BMI_online191008T204844.mat')
sum(data.selfHits)
27/2
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191002\NVI12\D7\BMI_online191002T151641.mat')
sum(data.selfHits)
clc
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191002\NVI13\D7\BMI_online191002T175532.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191002\NVI16\D7\BMI_online191002T204639.mat')
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191016\NVI12\D21\BMI_online191016T163233.mat')
sum(data.selfHits)
(34+12)/2
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191016\NVI13\D21\BMI_online191016T200634.mat')
sum(data.selfHits)
34/2
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191016\NVI16\D21\BMI_online191016T234334.mat')
sum(data.selfHits)
clear all
load('H:\holobmi_H\191025\NVI13\D25\BMI_online191025T181338.mat')
sum(data.selfHits)
clear all
[2 99]
y
n
[]
n
y
n
load('E:\holobmi_E\191025\NVI16\D25\holostim_seq191025T190913.mat')
uiopen('E:\holobmi_E\191025\NVI16\D25\im\holostim_seq\holostim_seq_191025T190815-098\holostim_seq_191025T190815-098_Cycle00001_VoltageRecording_001.csv',1)
%-- 10/25/2019 7:16 PM --%
load('G:\VivekNuria\Code\HoloBMI\analysis\nuria\checking_itworks.mat')
yy40 = [nanmean(nvi12_holo_40_total), nanmean(nvi12_rr_40_total), nanmean(nvi13_holo_40_total), nanmean(nvi13_rr_40_total), nanmean(nvi16_holo_40_total), nanmean(nvi16_rr_40_total)
)
yy40 = [nanmean(nvi12_holo_40_total), nanmean(nvi12_rr_40_total), nanmean(nvi13_holo_40_total), nanmean(nvi13_rr_40_total), nanmean(nvi16_holo_40_total), nanmean(nvi16_rr_40_total)]
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191018\NVI16\D23\BMI_online191018T202118.mat')
sum(data.selfHits)
nvi16_rr_40_total = [49]
yy40 = [nanmean(nvi12_holo_40_total), nanmean(nvi12_rr_40_total), nanmean(nvi13_holo_40_total), 0, nanmean(nvi16_holo_40_total), nanmean(nvi16_rr_40_total)]
err40 = [std(nvi12_holo_40_total)/sqrt(length(nvi12_holo_40_total)), std(nvi12_rr_40_total)/sqrt(length(nvi12_rr_40_total)), std(nvi13_holo_40_total)/sqrt(length(nvi13_holo_40_total)), 0 , std(nvi16_holo_40_total)/sqrt(length(nvi16_holo_40_total)), std(nvi16_rr_40_total)/sqrt(length(nvi16_rr_40_total))]
bar(xx,yy40)
errorbar(xx,yy40,err40)
hold on
errorbar(xx,yy40,err40)
bar(xx,yy40)
hold on
errorbar(xx,yy40,err40)
xticklabels({'1_H', '1_C', '2_H', '2_C', '3_H', '3_C'})
hline(20)
line(20)
hline(20)
yline(20)
line(20)
figure()
bar(xx,yy)
hold on
errorbar(xx,yy,err)
hline(20)
ylim([0,50])
xticklabels({'1_H', '1_C', '2_H', '2_C', '3_H', '3_C'})
title('15min')
title('40min')
nvi12hhits = nanmean(nv112_selfhits,1);
nvi12hhits = nanmean(nvi12_selfhits,1);
nvi13hhits = nanmean(nvi13_selfhits,1);
nvi16hhits = nanmean(nvi16_selfhits,1);
nvi12rrhits = nanmean(nvi12_rrselfhits,1);
nvi13rrhits = nanmean(nvi13_rrselfhits,1);
nvi16rrhits = nanmean(nvi16_rrselfhits,1);
nvi12hh = reshape(nvi12hhits,[1800,60]);
nvi13hh = reshape(nvi13hhits,[1800,60]);
nvi16hh = reshape(nvi16hhits,[1800,60]);
nvi12rr = reshape(nvi12rrhits,[1800,60]);
nvi13rr = reshape(nvi13rrhits,[1800,60]);
nvi16rr = reshape(nvi16rrhits,[1800,60]);
nvi12hpm = sum(nvi12hh, 1);
nvi13hpm = sum(nvi13hh, 1);
nvi16hpm = sum(nvi16hh, 1);
nvi12rrhpm = sum(nvi12rr, 1);
nvi13rrhpm = sum(nvi13rr, 1);
nvi16rrhpm = sum(nvi16rr, 1);
plot(nvi12hpm)
plot(smooth(nvi12hpm,2))
plot(smoothdata(nvi12hpm,2))
nvi12hpm(40:end)=[];
nvi13hpm(40:end)=[];
nvi16hpm(40:end)=[];
nvi12rrhpm(40:end)=[];
nvi13rrhpm(40:end)=[];
nvi16rrhpm(40:end)=[];
plot(smoothdata(nvi12hpm,2))
plot(smoothdata(nvi12rrhpm,2))
plot(smoothdata(nvi12hpm,3))
plot(smoothdata(nvi12hpm,2))
plot(smoothdata(nvi12hpm,3))
hold on
plot(smoothdata(nvi12hpm,5))
plot(smoothdata(nvi12hpm,'movmean',3))
plot(smoothdata(nvi12hpm,'movmean',2))
plot(smoothdata(nvi12hpm,'movmean',3))
plot(smoothdata(nvi12hpm,'movmean',1))
plot(smoothdata(nvi12hpm,'movmean',2))
hold on
plot(smoothdata(nvi12rrhpm,'movmean',2))
15/40
plot(smoothdata(nvi12rrhpm,'movmean',4))
plot(smoothdata(nvi12hpm,'movmean',5))
hold on
plot(smoothdata(nvi12rrhpm,'movmean',5))
title('animal 1')
legend({'H', 'C'})
xlabel('hpm (5minwindow)')
ylabel('hpm (5minwindow)')
xlabel('')
ylabel('hpm')
figure()
plot(smoothdata(nvi13hpm,'movmean',5))
hold on; plot(smoothdata(nvi13rrhpm,'movmean',5))
ylabel('hpm')
xlabel('Time (min)')
legend({'H', 'C'})
ylim([0.2, 1.4])
figure()
plot(smoothdata(nvi16hpm,'movmean',5))
hold on; plot(smoothdata(nvi16rrhpm,'movmean',5))
xlabel('Time (min)')
ylim([0.2, 1.4])
ylabel('hpm')
title('animal 3')
legend({'H', 'C'})
title('animal 2')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191015\NVI12\D20\BMI_online191015T143410.mat')
nvi12_40selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191015\NVI13\D20\BMI_online191015T175005.mat')
nvi13_40selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191015\NVI16\D20\BMI_online191015T204909.mat')
nvi16_40selfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191017\NVI13\D22\BMI_online191017T184758.mat')
nvi13_40selfhits = [nvi13_40selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191017\NVI16\D22\BMI_online191017T221038.mat')
nvi16_40selfhits = [nvi16_40selfhits; data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191018\NVI12\D23\BMI_online191018T152337.mat')
nvi12_40rrselfhits = [data.selfHits];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191018\NVI16\D23\BMI_online191018T202118.mat')
nvi16_40rrselfhits = [data.selfHits];
load('E:\holobmi_E\191025\NVI16\D25\BMI_online191025T195532.mat')
sum(data.holoHits)
nvi1240hhits = nanmean(nvi12_40selfhits,1);
nvi1340hhits = nanmean(nvi13_40selfhits,1);
nvi1640hhits = nanmean(nvi16_40selfhits,1);
nvi1240rrhhits = nanmean(nvi12_40rrselfhits,1);
nvi1640rrhhits = nanmean(nvi16_40rrselfhits,1);
nvi1240hh = reshape(nvi1240hhits,[1800,60]);
nvi1340hh = reshape(nvi1340hhits,[1800,60]);
nvi1640hh = reshape(nvi1640hhits,[1800,60]);
nvi1240rrhh = reshape(nvi1640rrhhits,[1800,60]);
nvi1640rrhh = reshape(nvi1640rrhhits,[1800,60]);
nvi1240hpm = sum(nvi1240hh, 1);
nvi1340hpm = sum(nvi1340hh, 1);
nvi1640hpm = sum(nvi1640hh, 1);
nvi1240rrhpm = sum(nvi1240rrhh, 1);
nvi1640rrhpm = sum(nvi1640rrhh, 1);
nvi1240hpm(40:end)=[];
nvi1340hpm(40:end)=[];
nvi1640hpm(40:end)=[];
nvi1240rrhpm(40:end)=[];
nvi1640rrhpm(40:end)=[];
close all
plot(smoothdata(nvi1240rrhpm,'movmean',5))
hold on; plot(smoothdata(nvi12rrhpm,'movmean',5))
plot(smoothdata(nvi1240rrhpm,'movmean',5))
hold on; plot(smoothdata(nvi1240rrhpm,'movmean',5))
plot(smoothdata(nvi1240rrhpm,'movmean',5))
hold on; plot(smoothdata(nvi1240rrhpm,'movmean',5))
plot(smoothdata(nvi1240hpm,'movmean',5))
hold on; plot(smoothdata(nvi1240rrhpm,'movmean',5))
xlabel('Time (min)')
ylabel('hpm')
title('animal 1')
legend({'H', 'C'})
figure()
plot(smoothdata(nvi1340hpm,'movmean',5))
title('animal 2')
xlabel('Time (min)')
ylabel('hpm')
figure()
plot(smoothdata(nvi1640hpm,'movmean',5))
hold on; plot(smoothdata(nvi1640rrhpm,'movmean',5))
plot(smoothdata(nvi1240rrhpm,'movmean',5))
nvi1240rrhh = reshape(nvi1240rrhhits,[1800,60]);
nvi1240rrhpm = sum(nvi1240rrhh, 1);
nvi1240rrhpm(40:end)=[];
figure()
plot(smoothdata(nvi1240hpm,'movmean',5))
hold on; plot(smoothdata(nvi1240rrhpm,'movmean',5))
xlabel('Time (min)')
ylabel('hpm')
title('animal 2')
title('animal 1')
legend({'H', 'C'})
ylim([0, 1.4])
figure()
plot(smoothdata(nvi1640hpm,'movmean',5))
hold on; plot(smoothdata(nvi1640rrhpm,'movmean',5))
ylabel('hpm')
xlabel('Time (min)')
title('animal 3')
legend({'H', 'C'})
sum(data.holoHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191018\NVI16\D23\BMI_online191018T202118.mat')
sum(data.selfHits)
sum(nvi16_40rrselfhits)
nvi1640hhits
sum(nvi1640hhits)
nvi1640hh = reshape(nvi1640hhits,[1800,60]);
sum(nvi1640hhits)
np.sum(nvi16_40selfhits)
sum(nvi16_40selfhits)
sum(sum(nvi16_40selfhits))
sum(sum(nvi16_40rrselfhits))
nvi1640hhits = nanmean(nvi16_40selfhits,1);
nvi1640rrhhits = nanmean(nvi16_40rrselfhits,1);
nvi1640hh = reshape(nvi1640hhits,[1800,60]);
nvi1640rrhh = reshape(nvi1640rrhhits,[1800,60]);
nvi1640hpm = sum(nvi1640hh, 1);
nvi1640rrhpm = sum(nvi1640rrhh, 1);
nvi1640hpm(40:end)=[];
nvi1640rrhpm(40:end)=[];
figure()
plot(smoothdata(nvi1640hpm,'movmean',5))
plot(smoothdata(nvi1640rrhpm,'movmean',5))
sum(nvi1640hhits)
sum(sum(nvi1640hh))
sum(sum(nvi1640hpm))
plot(nvi1640hhits)
plot(nvi1640rrhhits)
sum(nvi1640rrhhits)
sum(sum(nvi1640rrhpm))
sum(sum(nvi1240rrhpm))
close all
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI12\D5\BMI_online190930T152419.mat')
nvi12cursor = [data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI13\D5\BMI_online190930T180325.mat')
nvi13cursor = [data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI16\D5\BMI_online190930T202917.mat')
nvi16cursor = [data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI12\D6\BMI_online191001T152524.mat')
nvi12cursor = [nvi12cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI13\D6\BMI_online191001T180737.mat')
nvi13cursor = [nvi13cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_online191001T203900.mat')
nvi16cursor = [nvi16cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191003\NVI13\D8\BMI_online191003T175350.mat')
nvi13cursor = [nvi13cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI12\D10\BMI_online191005T155340.mat')
nvi12cursor = [nvi12cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI13\D10\BMI_online191005T183253.mat')
nvi13cursor = [nvi13cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI16\D10\BMI_online191005T213208.mat')
nvi16cursor = [nvi16cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI12\D12\BMI_online191007T153922.mat')
nvi12cursor = [nvi12cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI13\D12\BMI_online191007T202814.mat')
nvi13cursor = [nvi13cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI16\D12\BMI_online191007T233928.mat')
nvi16cursor = [nvi16cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI12\D14\BMI_online191009T153223.mat')
nvi12cursor = [nvi12cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI13\D14\BMI_online191009T192110.mat')
nvi13cursor = [nvi13cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI16\D14\BMI_online191009T220103.mat')
nvi16cursor = [nvi16cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI12\D16\BMI_online191011T151240.mat')
nvi12cursor = [nvi12cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI13\D16\BMI_online191011T175050.mat')
nvi13cursor = [nvi13cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI16\D16\BMI_online191011T205824.mat')
nvi16cursor = [nvi16cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191013\NVI13\D18\BMI_online191013T180942.mat')
nvi13cursor = [nvi13cursor; data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191013\NVI13\D18\BMI_target_info_20191013T164122.mat')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI12\D5\BMI_target_info_20190930T134821.mat')
nvi12T = [T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI13\D5\BMI_target_info_20190930T162928.mat')
nvi13T = [T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI16\D5\BMI_target_info_20190930T185737.mat')
nvi16T = [T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI12\D6\BMI_target_info_20191001T135333.mat')
nvi12T = [nvi12T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI13\D6\BMI_target_info_20191001T163732.mat')
nvi13T = [nvi13T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_target_info_20191001T191111.mat')
nvi16T = [nvi16T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191003\NVI13\D8\BMI_target_info_20191003T162208.mat')
nvi13T = [nvi13T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI12\D10\BMI_target_info_20191005T140501.mat')
nvi12T = [nvi12T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI13\D10\BMI_target_info_20191005T165116.mat')
nvi13T = [nvi13T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI16\D10\BMI_target_info_20191005T195923.mat')
nvi16T = [nvi16T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI12\D12\BMI_target_info_20191007T135419.mat')
nvi12T = [nvi12T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI13\D12\BMI_target_info_20191007T185256.mat')
nvi13T = [nvi13T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI16\D12\BMI_target_info_20191007T220344.mat')
nvi16T = [nvi16T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI12\D14\BMI_target_info_20191009T140400.mat')
nvi12T = [nvi12T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI13\D14\BMI_target_info_20191009T175422.mat')
nvi13T = [nvi13T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI16\D14\BMI_target_info_20191009T203233.mat')
nvi16T = [nvi16T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI12\D16\BMI_target_info_20191011T132219.mat')
nvi12T = [nvi12T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI13\D16\BMI_target_info_20191011T162138.mat')
nvi13T = [nvi13T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI16\D16\BMI_target_info_20191011T193113.mat')
nvi16T = [nvi16T, T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191013\NVI13\D18\BMI_target_info_20191013T164122.mat')
nvi13T = [nvi13T, T1];
nvi12T
length(find(nvicursor(1,:)>niv12T(1)))
length(find(nvi12cursor(1,:)>niv12T(1)))
length(find(nvi12cursor(1,:)>nvi12T(1)))
length(find(nvi12cursor(1,:)>nvi12T(1)))/length(nvi12cursor)
length(find(nvi12cursor(1,:)>nvi12T(1)))/length(nvi12cursor)*100
length(find(nvi12cursor(1,:)>nvi12T(1)))/length(find(nvi12cursor>0))*100
plt.plot(nvi12cursor)
plot(nvi12cursor)
plot(nvi12cursor(1,:))
length(find(nvi12cursor(1,:)>nvi12T(1)))/length(~isnan(nvi12cursor))*100
nvi12cursor(end-10:end)
length(find(nvi12cursor(1,:)>nvi12T(1)))/length(find(~isnan(nvi12cursor)))*100
length(find(~isnan(nvi12cursor))
length(find(~isnan(nvi12cursor)))
length(nvi12cursor)
length(find(nvi12cursor(1,:)>nvi12T(1)))/length(find(~isnan(nvi12cursor(1))))*100
length(find(nvi12cursor(1,:)>nvi12T(1)))/length(find(~isnan(nvi12cursor(1,:))))*100
length(find(~isnan(nvi12cursor(1))))
length(find(~isnan(nvi12cursor(1,:))))
length(find(nvi12cursor(1,:)>nvi12T(1)))/length(find(~isnan(nvi12cursor(1,:))))*100
for i=1:7
for i=1:7
nvi12hot = [];
for i=1:7
nvi12cold = [];
for i=1:7
nvi12hot[end] = length(find(nvi12cursor(i,:)>nvi12T(i)))/length(find(~isnan(nvi12cursor(i,:))))*100;
nvi12cold[end] = length(find(nvi12cursor(i,:)<-nvi12T(i)))/length(find(~isnan(nvi12cursor(i,:))))*100;
for i=1:6
nvi12hot(end) = length(find(nvi12cursor(i,:)>nvi12T(i)))/length(find(~isnan(nvi12cursor(i,:))))*100;
nvi12cold(end) = length(find(nvi12cursor(i,:)<-nvi12T(i)))/length(find(~isnan(nvi12cursor(i,:))))*100;
end
for i=1:7
nvi12hot(end+1) = length(find(nvi12cursor(i,:)>nvi12T(i)))/length(find(~isnan(nvi12cursor(i,:))))*100;
nvi12cold(end+1) = length(find(nvi12cursor(i,:)<-nvi12T(i)))/length(find(~isnan(nvi12cursor(i,:))))*100;
end
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI12\D6\BMI_online191001T152524.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
for i=1:8
nvi13hot(end+1) = length(find(nvi13cursor(i,:)>nvi13T(i)))/length(find(~isnan(nvi13cursor(i,:))))*100;
nvi13cold(end+1) = length(find(nvi13cursor(i,:)<-nvi13T(i)))/length(find(~isnan(nvi13cursor(i,:))))*100;
end
nvi12hot = [];
nvi12cold = [];
nvi13hot = [];
nvi13cold = [];
nvi16hot=[];
nvi16cold=[];
nvi12hotbb = [];
nvi12coldbb = [];
nvi13hotbb = [];
nvi13coldbb = [];
nvi16hotbb=[];
nvi16coldbb=[];

for i=1:6
    nvi12hot(end+1) = length(find(nvi12cursor(i,:)>=nvi12T(i)))/length(find(~isnan(nvi12cursor(i,:))))*100;
    nvi12cold(end+1) = length(find(nvi12cursor(i,:)<=-nvi12T(i)))/length(find(~isnan(nvi12cursor(i,:))))*100;
    nvi12hotbb(end+1) = length(find(nvi12baselinecursor(i,:)>=nvi12T(i)))/length(find(~isnan(nvi12baselinecursor(i,:))))*100;
    nvi12coldbb(end+1) = length(find(nvi12baselinecursor(i,:)<=-nvi12T(i)))/length(find(~isnan(nvi12baselinecursor(i,:))))*100;
end
for i=1:8
    nvi13hot(end+1) = length(find(nvi13cursor(i,:)>=nvi13T(i)))/length(find(~isnan(nvi13cursor(i,:))))*100;
    nvi13cold(end+1) = length(find(nvi13cursor(i,:)<=-nvi13T(i)))/length(find(~isnan(nvi13cursor(i,:))))*100;
    nvi13hotbb(end+1) = length(find(nvi13baselinecursor(i,:)>=nvi13T(i)))/length(find(~isnan(nvi13baselinecursor(i,:))))*100;
    nvi13coldbb(end+1) = length(find(nvi13baselinecursor(i,:)<=-nvi13T(i)))/length(find(~isnan(nvi13baselinecursor(i,:))))*100;
end
for i=1:6
    nvi16hot(end+1) = length(find(nvi16cursor(i,:)>=nvi16T(i)))/length(find(~isnan(nvi16cursor(i,:))))*100;
    nvi16cold(end+1) = length(find(nvi16cursor(i,:)<=-nvi16T(i)))/length(find(~isnan(nvi16cursor(i,:))))*100;
    nvi16hotbb(end+1) = length(find(nvi16baselinecursor(i,:)>=nvi16T(i)))/length(find(~isnan(nvi16baselinecursor(i,:))))*100;
    nvi16coldbb(end+1) = length(find(nvi16baselinecursor(i,:)<=-nvi16T(i)))/length(find(~isnan(nvi16baselinecursor(i,:))))*100;
end
nvi12gainh = (nvi12hot - nvi12hotbb)./nvi12hotbb*100;
nvi13gainh = (nvi13hot - nvi13hotbb)./nvi13hotbb*100;
nvi16gainh = (nvi16hot - nvi16hotbb)./nvi16hotbb*100;

nvi12gainc = (nvi12cold - nvi12coldbb)./nvi12hotbb*100;
nvi13gainc = (nvi13cold - nvi13coldbb)./nvi13hotbb*100;
nvi16gainc = (nvi16cold - nvi16coldbb)./nvi16hotbb*100;

gains = [nanmean(nvi12gainh), nanmean(nvi12gainc), nanmean(nvi13gainh), nanmean(nvi13gainc), nanmean(nvi16gainh), nanmean(nvi16gainc)]
errgains = [nanstd(nvi12gainh)/sqrt(length(nvi12gainh)), nanstd(nvi12gainc)/sqrt(length(nvi12gainc)), nanstd(nvi13gainh)/sqrt(length(nvi13gainh)), nanstd(nvi13gainc)/sqrt(length(nvi13gainc)), nanstd(nvi16gainh)/sqrt(length(nvi16gainh)), nanstd(nvi16gainc)/sqrt(length(nvi16gainc))]

nvi12hotrr = [];
nvi12coldrr = [];
nvi13hotrr = [];
nvi13coldrr = [];
nvi16hotrr=[];
nvi16coldrr=[];
nvi12hotrrbb = [];
nvi12coldrrbb = [];
nvi13hotrrbb = [];
nvi13coldrrbb = [];
nvi16hotrrbb=[];
nvi16coldrrbb=[];

for i=1:3
    nvi12hotrr(end+1) = length(find(nvi12cursorrr(i,:)>=nvi12Trr(i)))/length(find(~isnan(nvi12cursorrr(i,:))))*100;
    nvi12coldrr(end+1) = length(find(nvi12cursorrr(i,:)<=-nvi12Trr(i)))/length(find(~isnan(nvi12cursorrr(i,:))))*100;
    nvi12hotrrbb(end+1) = length(find(nvi12baselinecursorrr(i,:)>=nvi12Trr(i)))/length(find(~isnan(nvi12baselinecursorrr(i,:))))*100;
    nvi12coldrrbb(end+1) = length(find(nvi12baselinecursorrr(i,:)<=-nvi12Trr(i)))/length(find(~isnan(nvi12baselinecursorrr(i,:))))*100;
end
for i=1:3
    nvi13hotrr(end+1) = length(find(nvi13cursorrr(i,:)>=nvi13Trr(i)))/length(find(~isnan(nvi13cursorrr(i,:))))*100;
    nvi13coldrr(end+1) = length(find(nvi13cursorrr(i,:)<=-nvi13Trr(i)))/length(find(~isnan(nvi13cursorrr(i,:))))*100;
    nvi13hotrrbb(end+1) = length(find(nvi13baselinecursorrr(i,:)>=nvi13Trr(i)))/length(find(~isnan(nvi13baselinecursorrr(i,:))))*100;
    nvi13coldrrbb(end+1) = length(find(nvi13baselinecursorrr(i,:)<=-nvi13Trr(i)))/length(find(~isnan(nvi13baselinecursorrr(i,:))))*100;
end
for i=1:3
    nvi16hotrr(end+1) = length(find(nvi16cursorrr(i,:)>=nvi16Trr(i)))/length(find(~isnan(nvi16cursorrr(i,:))))*100;
    nvi16coldrr(end+1) = length(find(nvi16cursorrr(i,:)<=-nvi16Trr(i)))/length(find(~isnan(nvi16cursorrr(i,:))))*100;
    nvi16hotrrbb(end+1) = length(find(nvi16baselinecursorrr(i,:)>=nvi16Trr(i)))/length(find(~isnan(nvi16baselinecursorrr(i,:))))*100;
    nvi16coldrrbb(end+1) = length(find(nvi16baselinecursorrr(i,:)<=-nvi16Trr(i)))/length(find(~isnan(nvi16baselinecursorrr(i,:))))*100;
end
nvi12gainrrh = (nvi12hotrr - nvi12hotrrbb)./nvi12hotrrbb*100;
nvi13gainrrh = (nvi13hotrr - nvi13hotrrbb)./nvi13hotrrbb*100;
nvi16gainrrh = (nvi16hotrr - nvi16hotrrbb)./nvi16hotrrbb*100;

nvi12gainrrc = (nvi12coldrr - nvi12coldrrbb)./nvi12hotrrbb*100;
nvi13gainrrc = (nvi13coldrr - nvi13coldrrbb)./nvi13hotrrbb*100;
nvi16gainrrc = (nvi16coldrr - nvi16coldrrbb)./nvi16hotrrbb*100;

gainsrr = [nanmean(nvi12gainrrh), nanmean(nvi12gainrrc), nanmean(nvi13gainrrh), nanmean(nvi13gainrrc), nanmean(nvi16gainrrh), nanmean(nvi16gainrrc)]
errgainsrr = [nanstd(nvi12gainrrh)/sqrt(length(nvi12gainrrh)), nanstd(nvi12gainrrc)/sqrt(length(nvi12gainrrc)), nanstd(nvi13gainrrh)/sqrt(length(nvi13gainrrh)), nanstd(nvi13gainrrc)/sqrt(length(nvi13gainrrc)), nanstd(nvi16gainrrh)/sqrt(length(nvi16gainrrh)), nanstd(nvi16gainrrc)/sqrt(length(nvi16gainrrc))]



xx = [1,2, 4,5, 7,8];
yyhot = [nanmean(nvi12hot), nanmean(nvi12cold), nanmean(nvi13hot), nanmean(nvi13cold), nanmean(nvi16hot), nanmean(nvi16cold)]
errhot = [std(nvi12hot)/sqrt(length(nvi12hot)), std(nvi12cold)/sqrt(length(nvi12cold)), std(nvi13hot)/sqrt(length(nvi13hot)), std(nvi13cold)/sqrt(length(nvi13cold)), std(nvi16hot)/sqrt(length(nvi16hot)), std(nvi16cold)/sqrt(length(nvi16cold))]
figure()
plot(xx,yyhot)
bar(xx,yyhot)
errorbar(xx,yyhot,errhot)
bar(xx,yyhot)
hold on
errorbar(xx,yyhot,errhot)
xticklabels({'1_H', '1_C', '2_H', '2_C', '3_H', '3_C'})
ylabel('%occupancy')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191004\NVI12\D9\BMI_online191004T153616.mat')
nvi12cursorE3 = [data.cursor];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191004\NVI12\D9\BMI_target_info_20191004T145314.mat')
nvi12TE3 = [T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191004\NVI13\D9\BMI_online191004T183145.mat')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191004\NVI13\D9\BMI_target_info_20191004T174606.mat')
nvi13cursorE3 = [data.cursor];
nvi13TE3 = [T1];
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191004\NVI16\D9\BMI_online191004T212633.mat')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191004\NVI16\D9\BMI_target_info_20191004T204004.mat')
nvi16cursorE3 = [data.cursor];
nvi16TE3 = [T1];
nvi16hot = length(find(nvi16cursorE3>nvi16TE3))/length(find(~isnan(nvi16cursorE3)))*100;
nvi13hot = length(find(nvi13cursorE3>nvi13TE3))/length(find(~isnan(nvi13cursorE3)))*100;
nvi12hot = length(find(nvi12cursorE3>nvi12TE3))/length(find(~isnan(nvi12cursorE3)))*100;
nvi16cold = length(find(nvi16cursorE3<-nvi16TE3))/length(find(~isnan(nvi16cursorE3)))*100;
nvi13cold = length(find(nvi13cursorE3<-nvi13TE3))/length(find(~isnan(nvi13cursorE3)))*100;
nvi12cold = length(find(nvi12cursorE3<-nvi12TE3))/length(find(~isnan(nvi12cursorE3)))*100;
nvi12hot
nvi12cold
yyhotE3 = [nanmean(nvi12hot), nanmean(nvi12cold), nanmean(nvi13hot), nanmean(nvi13cold), nanmean(nvi16hot), nanmean(nvi16cold)]
figure()
bar(xx,yyhotE3)
hold on
ylabel('%occupancy')
xticklabels({'1_H', '1_C', '2_H', '2_C', '3_H', '3_C'})
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191004\NVI16\D9\BMI_online191004T221010.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191004\NVI12\D9\BMI_online191004T153616.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI12\D6\BMI_online191001T152524.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_online191001T203900.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI13\D5\BMI_online190930T180325.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI16\D10\BMI_online191005T213208.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI16\D10\BMI_online191005T221448.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI13\D5\BMI_online190930T180325.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI12\D6\BMI_online191001T152524.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_online191001T203900.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_online191001T212200.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI12\D10\BMI_online191005T155340.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI16\D10\BMI_online191005T213208.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191005\NVI16\D10\BMI_online191005T221448.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI12\D12\BMI_online191007T153922.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI12\D12\BMI_online191007T155538.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI13\D12\BMI_online191007T202814.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI16\D12\BMI_online191007T233928.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI12\D14\BMI_online191009T153223.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI12\D14\BMI_online191009T161633.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
figure
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI13\D14\BMI_online191009T192110.mat')
figure
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI16\D14\BMI_online191009T220103.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI12\D16\BMI_online191011T151240.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI12\D16\BMI_target_info_20191011T132219.mat')
T1
figure
kk = nvi12_selfhits(6,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi12_selfhits(1,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi12_selfhits(2,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi12_selfhits(3,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi12_selfhits(4,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi12_selfhits(5,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi16_selfhits(1,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi16_selfhits(2,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi16_selfhits(3,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi16_selfhits(4,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi16_selfhits(5,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi16_selfhits(6,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
kk = nvi16_selfhits(2,:);
kk2 = reshape(kk,[1800,60]);
kk3 = sum(kk2, 1);
kk3(41:end) = [];
plot(smoothdata(kk3,'movmean',5))
plot(smoothdata(kk3,'movmean',8))
plot(smoothdata(kk3,'movmean',2))
plot(smoothdata(kk3,'movmean',5))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191001\NVI16\D6\BMI_online191001T212200.mat')
kk = data.selfHits;
kk2 = reshape(kk,[1800,60]);
kk3 = [kk3, sum(kk2, 1)];
kk3(81:end) = [];
plot(smoothdata(kk3,'movmean',5))
for i=1:7
nvi12hot(end+1) = length(find(nvi12cursor(i,:)>nvi12T(i)*8/10))/length(find(~isnan(nvi12cursor(i,:))))*100;
nvi12cold(end+1) = length(find(nvi12cursor(i,:)<-nvi12T(i)*8/10))/length(find(~isnan(nvi12cursor(i,:))))*100;
end
nvi12cold = [];
nvi12hot = [];
for i=1:6
nvi12hot(end+1) = length(find(nvi12cursor(i,:)>nvi12T(i)*8/10))/length(find(~isnan(nvi12cursor(i,:))))*100;
nvi12cold(end+1) = length(find(nvi12cursor(i,:)<-nvi12T(i)*8/10))/length(find(~isnan(nvi12cursor(i,:))))*100;
end
nvi12hot
nvi12cold
nvi13cold = [];
nvi13hot = [];
nvi16hot=[]
nvi16cold=[]
for i=1:8
nvi13hot(end+1) = length(find(nvi13cursor(i,:)>nvi13T(i)*8/10))/length(find(~isnan(nvi13cursor(i,:))))*100;
nvi13cold(end+1) = length(find(nvi13cursor(i,:)<-nvi13T(i)*8/10))/length(find(~isnan(nvi13cursor(i,:))))*100;
end
for i=1:6
nvi16hot(end+1) = length(find(nvi16cursor(i,:)>nvi16T(i)*8/10))/length(find(~isnan(nvi16cursor(i,:))))*100;
nvi16cold(end+1) = length(find(nvi16cursor(i,:)<-nvi16T(i)*8/10))/length(find(~isnan(nvi16cursor(i,:))))*100;
end
yyhotE3 = [nanmean(nvi12hot), nanmean(nvi12cold), nanmean(nvi13hot), nanmean(nvi13cold), nanmean(nvi16hot), nanmean(nvi16cold)]
yyhot = [nanmean(nvi12hot), nanmean(nvi12cold), nanmean(nvi13hot), nanmean(nvi13cold), nanmean(nvi16hot), nanmean(nvi16cold)]
bar(xx,yyhot)
errhot = [std(nvi12hot)/sqrt(length(nvi12hot)), std(nvi12cold)/sqrt(length(nvi12cold)), std(nvi13hot)/sqrt(length(nvi13hot)), std(nvi13cold)/sqrt(length(nvi13cold)), std(nvi16hot)/sqrt(length(nvi16hot)), std(nvi16cold)/sqrt(length(nvi16cold))]
hold on
errorbar(xx,yyhot,errhot)
ylabel('%occupancy')
title('20% to T1')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191024\NVI12\D24\BMI_online191024T171438.mat')
sum(sum(nvi1240rrhpm))
sum(data.selfHits)
clear all
%-- 10/26/2019 11:20 AM --%
load('G:\VivekNuria\Code\HoloBMI\analysis\nuria\checking_itworks_ruiquestions.mat')
plot(smoothdata(nvi12hpm,'movmean',5))
hpm = [nvi12hpm; nvi13hpm; nvi16hpm];
rrhpm = [nvi12rrhpm; nvi13rrhpm; nvi16rrhpm];
plot(nanmean(hpm,1))
hold on
plot(nanmean(rrhpm,1))
plot(smoothdata(nanmean(hpm,1), 'movmean',5)
plot(smoothdata(nanmean(hpm,1), 'movmean',5))
hold on
plot(smoothdata(nanmean(rrhpm,1), 'movmean',5))
figure
hpmall = [nvi12hh; nvi13hh;nvi16hh];
hpmrrall = [nvi12rrhh; nvi13rrhh;nvi16rrhh];
hpmrrall = [nvi12rr; nvi13rr;nvi16rr];
plot(smoothdata(nanmean(hpmall,1), 'movmean',5))
hpmall(40:end)=[]
hpmall = [nvi12hh; nvi13hh;nvi16hh];
hpmall = [nvi12_selfhits; nvi13_selfhits; nvi16_selfhits];
hpmalll = reshape(np.nanmean(hpmall,1),[1800,60]);
hpmalll = reshape(nanmean(hpmall,1),[1800,60]);
kk = nanmean(hpmall,1);
1800*60
hpmalll = reshape(nanmean(hpmall,1),[1800,60]);
hpmalll = reshape(kk,[1800,60]);
hpmalll = reshape(nanmean(hpmall,1),[1800,60]);
hpm_tog = sum(hpmalll, 1);
hpm_tog(40:end)=[];
hpmallrr = [nvi12_rrselfhits; nvi13_rrselfhits; nvi16_rrselfhits];
hpmarrll = reshape(nanmean(hpmallrr,1),[1800,60]);
hpm_rrtog = sum(hpmarrlll, 1);
hpm_rrtog = sum(hpmarrll, 1);
hpm_rrtog(40:end)=[];
plt.plot(hpm_tog)
plot(hpm_tog)
plot(smoothdata(hpm_tog, 'movmean',5))
hold on
plot(smoothdata(hpm_rrtog, 'movmean',5))
xlabel('Time (min)')
ylabel('hpm')
legend({'H', 'C'})
yy
15total = [nanmean(nvi12_15_total), nanmean(nvi13_15_total), nanmean(nvi16_15_total)]
total = [nanmean(nvi12_15_total), nanmean(nvi13_15_total), nanmean(nvi16_15_total)]
holo15 = [nvi12_holo_15_total, nvi13_holo_15_total, nvi16_holo_15_total];
holo15av = [nanmean(nvi12_holo_15_total), nanmean(nvi13_holo_15_total), nanmean(nvi16_holo_15_total)];
rr15 = [nvi12_rr_15_total, nvi13_rr_15_total, nvi16_rr_15_total];
rr15av = [nanmean(nvi12_rr_15_total), nanmean(nvi13_rr_15_total), nanmean(nvi16_rr_15_total)];
yy15 = [nanmean(holo15), nanmean(rr15)]
yy15av = [nanmean(holo15av), nanmean(rr15av)]
bar([1,2], yy15av)
rr15 = [nanstd(holo15)/sqrt(length(holo15)), nanstd(rr15)/sqrt(length(rr15))]
rr15av = [nanstd(holo15av)/sqrt(length(holo15av)), nanstd(rr15av)/sqrt(length(rr15av))]
hold on
errorbar([1,2],yy15av,err15av)
errorbar([1,2],yy15av,rr15av)
rr15 = [nvi12_rr_15_total, nvi13_rr_15_total, nvi16_rr_15_total];
rr15av = [nanmean(nvi12_rr_15_total), nanmean(nvi13_rr_15_total), nanmean(nvi16_rr_15_total)];
error15av = [nanmean(nvi12_rr_15_total), nanmean(nvi13_rr_15_total), nanmean(nvi16_rr_15_total)];
error15av = [nanstd(holo15av)/sqrt(length(holo15av)), nanstd(rr15av)/sqrt(length(rr15av))]
ttest(holo15av, rr15av)
[h,p] = ttest(holo15av, rr15av)
holo15av
rr15av
xticklabels({'H', 'C'})
ylabel('Total hits')
bar([1,2], yy15av)
plot([1,2], [holo15av(1), rr15av(1))
plot([1,2], [holo15av(1), rr15av(1)])
bar([1,2], yy15av)
hold on
plot([1,2], [holo15av(1), rr15av(1)])
plot([1,2], [holo15av(2), rr15av(2)])
plot([1,2], [holo15av(1), rr15av(1)],'k')
plot([1,2], [holo15av(2), rr15av(2)],'k')
plot([1,2], [holo15av(3), rr15av(3)],'k')
err15av
error15av
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI12\D25\BMI_online191025T154101.mat')
prefdir
close all
clear all
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI12\D25\BMI_online191025T154101.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI13\D25\BMI_online191025T181338.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI16\D25\BMI_online191025T214847.mat')
bar([-0.5:0.01:0.9], histc(data.cursor, [-0.5:0.01:0.9]))
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI12\D16\BMI_online191011T151240.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI13\D16\BMI_online191011T175050.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI16\D16\BMI_online191011T205824.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
figure()
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI12\D14\BMI_online191009T153223.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI13\D14\BMI_online191009T192110.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
figure()
sum(data.selfHits)
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI16\D14\BMI_online191009T220103.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI12\D12\BMI_online191007T153922.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI13\D12\BMI_online191007T202814.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
sum(data.selfHits)
title('During BMI')
xlabel('Cursor dist')
ylabel('Counts')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191007\NVI16\D12\BMI_online191007T233928.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI13\D5\BMI_online190930T180325.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
title('During BMI')
xlabel('Cursor dist')
ylabel('Counts')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
figure
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI12\D25\BMI_online191025T154101.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI13\D25\BMI_online191025T181338.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI12\D16\BMI_online191011T151240.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI13\D16\BMI_online191011T175050.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191011\NVI16\D16\BMI_online191011T205824.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI12\D14\BMI_online191009T153223.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191009\NVI13\D14\BMI_online191009T192110.mat')
bar([-0.9:0.01:0.9], histc(data.cursor, [-0.9:0.01:0.9]))
sum(data.selfHits)
title('During BMI')
xlabel('Cursor dist')
ylabel('Counts')
40/15
clear all
load('I:\Vivek\Imaging_data\191025\NVI16\D25\BaselineOnline191025T191404.mat')
load('I:\Vivek\Imaging_data\191025\NVI16\D25\BMI_online191025T214847.mat')
clear all
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI16\D25\BMI_target_info_20191025T193334.mat')
num_E1 = length(E1_base);
num_E2 = length(E2_base);
num_neurons = num_E1 + num_E2;
E_id = [1*ones(num_E1, 1); 2*ones(num_E2, 1)];
E1_sel = E_id==1;
E2_sel = E_id==2;
E1_proj = zeros(num_neurons, 1);
E1_proj(E1_sel) = 1;
E1_norm = sum(E1_sel); %can replace with vector norm.
disp('E1 proj');
E1_proj = E1_proj/E1_norm;
E1_proj
E2_proj = zeros(num_neurons, 1);
E2_proj(E2_sel) = 1;
E2_norm = sum(E2_sel);
disp('E2 proj')
E2_proj = E2_proj/E2_norm;
E2_proj
disp('decoder:')
decoder = E2_proj - E1_proj
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI16\D25\BaselineOnline191025T191404.mat')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\191025\NVI16\D25\target_calibration_ALL_20191025T193334.mat')
bar([-0.9:0.01:0.9], histc(cursor_obs, [-0.9:0.01:0.9]))
clear all
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI12\D5\BMI_online190930T152419.mat')
load('Z:\Vivek\VivekNuria\HoloBMI_092319\190930\NVI12\D5\target_calibration_ALL_20190930T134821.mat')
T1