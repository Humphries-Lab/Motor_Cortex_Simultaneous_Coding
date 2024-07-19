% add the path to matnwb and generate the core classes
addpath('C:\Users\controlmotor\Desktop\Andrea\codes_from_papers\matnwb-master')
addpath('matnwb');
 
% Reminder: YOU DO NOT NORMALLY NEED TO CALL THIS FUNCTION. Only attempt this method if you
% encounter read errors.
generateCore(util.getSchemaVersion('sub-M_ses-CO-20140203_behavior+ecephys.nwb'));

% ignorecache informs the `nwbRead` call to not generate files by default. Since we have already
% done this, we can skip this automated step when reading. If you are reading the file before
% generating, you can omit this argument flag.
nwb = nwbRead('sub-M_ses-CO-20140203_behavior+ecephys.nwb', 'ignorecache')