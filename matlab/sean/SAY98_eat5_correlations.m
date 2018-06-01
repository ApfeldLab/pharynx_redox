%% Load Data
m_410 = csvread('../data/SAY98_eat5_2do_05_24/eat5_2do_05_24_410_measure.csv', 1, 1);
m_470 = csvread('../data/SAY98_eat5_2do_05_24/eat5_2do_05_24_470_measure.csv', 1, 1);

%% Register
[fdObjs, data] = ChannelRegister(ssquare(m_410), ssquare(m_470));

resample_size = size(data,1);

%%
data.m410 = struct('all', data.m410, 'regionMeans', regionMeans(data.m410));
data.m470 = struct('all', data.m470, 'regionMeans', regionMeans(data.m470));

data.R.all = data.m410.all ./ data.m470.all;
data.R.regionMeans = regionMeans(data.R.all);

data.OxD.all = ja_oxd(data.R.all);
data.OxD.regionMeans = regionMeans(data.OxD.all);

data.E.all = ja_E(data.OxD.all);
data.E.regionMeans = regionMeans(data.E.all);

data.E_minus_medial.all = data.E.all - (data.E.regionMeans.medial_axis * ones(1,1000)).';
data.E_minus_medial.regionMeans = regionMeans(data.E_minus_medial.all);

data.E_minus_pm3.all = data.E.all - (data.E.regionMeans.pm3 * ones(1,1000)).';
data.E_minus_pm3.regionMeans = regionMeans(data.E_minus_pm3.all);

data.E_minus_pm5.all = data.E.all - (data.E.regionMeans.pm5 * ones(1,1000)).';
data.E_minus_pm5.regionMeans = regionMeans(data.E_minus_pm5.all);

data.E_minus_pm7.all = data.E.all - (data.E.regionMeans.pm7 * ones(1,1000)).';
data.E_minus_pm7.regionMeans = regionMeans(data.E_minus_pm7.all);

%%
data_col_labels ={ 'i410','i470','R','OxD','E (mV)','E-Emedialaxis (mV)','E-Epm3 (mV)','E-Epm5 (mV)','E-Epm7 (mV)'};
positional_data = cat(3,int_410(:,:),int_470(:,:),R(:,:), ...
    OxD(:,:),E(:,:),...
    E(:,:)-(ja_E(ja_oxd(R_region(:,1)* ones(1,100))))' ,...
    E(:,:)-(ja_E(ja_oxd(R_region(:,5)* ones(1,100))))' ,...
    E(:,:)-(ja_E(ja_oxd(R_region(:,4)* ones(1,100))))',...
    E(:,:)-(ja_E(ja_oxd(R_region(:,2)* ones(1,100))))');
