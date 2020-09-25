
names = {'autodiff', ...
         'internal', ...
         'multiscale', ...
         'visualization', ...
         'model-io', ...
         'solvers'};
names = cellfun(@(x) fullfile(ROOTDIR, '..', ['mrst-', x]), names, ...
                    'UniformOutput', false);

mrstPath('addroot', names{:});
mrstPath('add', 'co2lab', fullfile(ROOTDIR, '..', 'mrst-co2lab', 'co2lab'));

mrstPath('add', 'FizedPointSolversRichards', '/Users/davide/Documents/GitHub/MRST_UiB/mrst-davideGitHub');
mrstModule add
clear namesmrstPath reregister matlab_bgl ...
    '/Users/davide/Documents/GitHub/MRST_UiB'
