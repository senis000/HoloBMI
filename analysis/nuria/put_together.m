function put_together()

animal(1).name = 'NVI12';
animal(2).name = 'NVI13';
animal(3).name = 'NVI16';
animal(4).name = 'NVI17';
animal(5).name = 'NVI20';
animal(6).name = 'NVI22';
%just to bring the values of the days when we did the experiments
animal(1).days.holobmi = [190930, 191001, 191005, 191007, 191009, 191011, 191013, 191025]; 
animal(2).days.holobmi = [190930, 191001, 191003, 191005, 191007, 191009, 191011, 191013, 191025]; 
animal(3).days.holobmi = [190930, 191001, 191005, 191007, 191009, 191011, 191025];

animal(1).days.E3 = [191004, 191026, 191028];
animal(2).days.E3 = [191004, 191026, 191028];
animal(3).days.E3 = [191004, 191026, 191028];

animal(1).days.rr = [191006, 191010, 191014];
animal(2).days.rr = [191006, 191010, 191014];
animal(3).days.rr = [191006, 191010, 191014, 191024];

animal(1).days.E2holo = [191008, 191012, 191027];
animal(2).days.E2holo = [191008, 191012, 191027];
animal(3).days.E2holo = [191008, 191012, 191027];

animal(1).days.frr = [191030, 191101, 191103];
animal(2).days.frr = [191030, 191101, 191103];
animal(3).days.frr = [191030, 191101, 191103];

animal(1).days.fhbmi = [191031, 191102, 191104];
animal(2).days.fhbmi = [191031, 191102, 191104];
animal(3).days.fhbmi = [191031, 191102, 191104];


animal(4).days.holobmi = [191106, 191108, 191110, 191112, 191114, 191116, 191120, 191122];
animal(5).days.holobmi = [191106, 191108, 191110, 191112, 191114, 191116, 191118, 191120];
animal(6).days.holobmi = [191106, 191108, 191110, 191114, 191116, 191118, 191120, 191122];

animal(4).days.E3 = [191109, 191115, 191121];
animal(5).days.E3 = [191109, 191115, 191121];
animal(6).days.E3 = [191109, 191115, 191121];

animal(4).days.rr = [191111, 191117, 121123];
animal(5).days.rr = [191105, 191117, 121122];
animal(6).days.rr = [191105, 191111, 191117];

animal(4).days.E2holo = [191107, 191113, 191119];
animal(5).days.E2holo = [191107, 191113, 191119];
animal(6).days.E2holo = [191107, 191113, 191119];

animal(4).days.frr = [191124, 191126, 191128];
animal(5).days.frr = [191124, 191126, 191128];
animal(6).days.frr = [191124, 191126, 191128];

animal(4).days.fhbmi = [191125, 191127, 191129];
animal(5).days.fhbmi = [191125, 191127, 191130];
animal(6).days.fhbmi = [191125, 191127, 191129];

%preparing the structure
dirpath = 'C:/Data/HoloBMI/experiments/';
for aa=1:length(animal)
    typeDays = fieldnames(animal(aa).days);
    for tt=1:length(typeDays)
        daystoGet = animal(aa).days.(typeDays{tt});
        for dd = daystoGet
            %this is all to find the bmi_online files
            folderPath = fullfile(dirpath,num2str(dd),animal(aa).name);
            filesPath = dir(folderPath);
            ibmi=1;
            itar=1;
            for ff=3:length(filesPath)
                if strcmp(filesPath(ff).name(1:5),'BMI_o')
                    bmiFiles{ibmi} = filesPath(ff).name;
                    ibmi = ibmi+1;
                end
                if strcmp(filesPath(ff).name(1:5),'targe')
                    targetFiles{itar} = filesPath(ff).name;
                    itar = itar+1;
                end
            end
            %if there is less than 2 bmi_online something is wrong
            if length(bmiFiles)<2
                disp('STHAAAAP Something is wrong with this experiment')
                disp(folderPath)
                return
            else
                fileBMIOpen = fullfile(folderPath, bmiFiles{2});
                fileTarOpen = fullfile(folderPath, targetFiles{end});
                load(fileBMIOpen, 'data');
                load(fileTarOpen, 'cursor_obs', 'hits', 'hits_valid');
            end
        end
    end
end



