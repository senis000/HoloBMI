% animals and dates to adjust

% animal= 'NVI12';
% date = '191007';
% BMI_online191007T153922 -> data_1
% BMI_online191007T155538 -> data_2
folder = 'C:/Users/Nuria/Documents/DATA/holoBMI/raw_data/191007/NVI12/D12/';
data_2 = cut_mat_files(data_2, 75600 - 63000); % to fill to 75600
merge_mat_files (data_1, data_2, bData_1, folder);


% animal= 'NVI12';
% date = '191025';
% BMI_online191025T143753 -> data_1
% BMI_online191025T145556 -> data_2
folder = 'C:/Users/Nuria/Documents/DATA/holoBMI/raw_data/191025/NVI12/D25/';
merge_mat_files (data_1, data_2, bData_1, folder);


% animal= 'NVI16';
% date = '191025';
% BMI_online191025T195532 -> data_1
% BMI_online191025T210459 -> data_2
folder = 'C:/Users/Nuria/Documents/DATA/holoBMI/raw_data/191025/NVI16/D25/';
merge_mat_files (data_1, data_2, bData_1, folder);


% animal= 'NVI16';
% date = '191026';
% BMI_online191026T195005 -> data_1
% BMI_online191026T204727 -> data_2
folder = 'C:/Users/Nuria/Documents/DATA/holoBMI/raw_data/191026/NVI16/D26/';
merge_mat_files (data_1, data_2, bData_1, folder);



% animal= 'NVI16';
% date = '191031';
% BMI_online191031T200328 -> data_1
% BMI_online191031T202659 -> data_2
% BMI_online191031T205744 -> data_3
% BMI_online191031T211608 -> data_4
folder = 'C:/Users/Nuria/Documents/DATA/holoBMI/raw_data/191031/NVI16/D31/';
merge_mat_files (data_1, data_2, bData_1, folder);
merge_mat_files (data_3, data_4, bData_3, folder);



