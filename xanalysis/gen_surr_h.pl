#!/usr/bin/perl

for (`cat x2000parameter.defs`) {
    next if /^[cd\#\*]/i;
    $in_param = 1 if /parameter/i;
    if ($in_param && /[^a-zA-Z_0-9]MAX_NUM_CHAN\s*=\s*(\d+)/) {
        $max_num_chan = $1;
    }
    if ($in_param && /[^a-zA-Z_0-9]MAX_NUM_CODES\s*=\s*(\d+)/) {
        $max_num_codes = $1;
    }
    if ($in_param && /[^a-zA-Z_0-9]iTOTAL_SHIFTS\s*=\s*(\d+)/) {
        $itotal_shifts = $1;
    }
    $in_param = 1 if /\)/i;
}

print <<EOT;
#define MAX_NUM_CHAN $max_num_chan
#define MAX_NUM_CODES $max_num_codes
#define iTOTAL_SHIFTS $itotal_shifts

void ref_surrogate (double *spiketimes, int *ital, int *idp, void *ssp, int *mne, int *spikecount_array);

void
tar_surrogate (void *spk_tm_ptr, int *mne, int *ital, int *ids,
               int *excluded, unsigned long *stp, int *surr_ital,
               double *E_begin, double *E_end, int *ecnt);

EOT
