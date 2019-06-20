
Step 1)
Run Associate_tracks_20151221.m to associate GFP and RFP tracks. (used to be 0909)
Run Associate_tracks_um_20150706.m to associate GFP and RFP tracks, if um unit was used for Imaris analysis.
(When you need to apply updated code to a pre-existing data, use Associate_update_from_master2ref_format.m)

Step 2)
Run Clean_associated_tracks2015_03_09.m 
or Clean_associated_tracks_batch.m for automatic cleaning
or Clean_associated_tracks_um_batch.m for automatic cleaning where um unit was used for Imaris analysis.
Run Manually_pick_asso_ref_cor20150312.m for manual cleaning.
Step 3)
Run Plot_stats20151006.m (was 20150910)

Useful short cuts:
shift+apple+d will open the file where cursor is on.
ctrl + C will abort the process. 