%% denoising
home = "/n/cohen_lab/Lab/Labmembers/Michael Xie/104438_Tread-1speed_Dilas-8V_488-OD0.0";
mov_in = "movReg.bin";
output = fullfile(home,"output");
detr_spacing = 500;
row_blocks = 6;
col_blocks = 2;

run_command = sprintf("source setup.sh\n sbatch denoise.run ""%s"" ""%s"" ""%s"" %d %d %d",...
    home, mov_in, output, detr_spacing, row_blocks, col_blocks);

system(run_command)