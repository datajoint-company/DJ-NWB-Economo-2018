import re
import os
import sys
from datetime import datetime

import numpy as np
import scipy.io as sio
import datajoint as dj
import tqdm

from . import (reference, utilities, acquisition, analysis)

schema = dj.schema(dj.config['custom'].get('database.prefix', '') + 'extracellular')


@schema
class ProbeInsertion(dj.Manual):
    definition = """ # Description of probe insertion details during extracellular recording
    -> acquisition.Session
    -> reference.Probe
    -> reference.BrainLocation
    insertion_depth: decimal(6,2)  #  (um)
    """


@schema
class UnitSpikeTimes(dj.Manual):
    definition = """ 
    -> ProbeInsertion
    unit_id : smallint
    ---
    -> reference.Probe.Channel
    unit_cell_type: enum('PTlower', 'PTupper')  #  depending on the animal (which cell-type being tagged)
    unit_quality="": varchar(32)  # 
    unit_depth: float  # (um)
    spike_times: longblob  # (s) time of each spike, with respect to the start of session 
    """


@schema
class TrialSegmentedUnitSpikeTimes(dj.Computed):
    definition = """
    -> UnitSpikeTimes
    -> acquisition.TrialSet.Trial
    -> analysis.TrialSegmentationSetting
    ---
    segmented_spike_times: longblob
    """

    def make(self, key):
        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')
        # get event time
        try:
            event_time_point = analysis.get_event_time(event_name, key)
        except analysis.EventChoiceError as e:
            print(f'Trial segmentation error - Msg: {str(e)}')
            return

        pre_stim_dur = float(pre_stim_dur)
        post_stim_dur = float(post_stim_dur)

        # get raw & segment
        spike_times = (UnitSpikeTimes & key).fetch1('spike_times')
        key['segmented_spike_times'] = spike_times[np.logical_and(
            (spike_times >= (event_time_point - pre_stim_dur)),
            (spike_times <= (event_time_point + post_stim_dur)))] - event_time_point

        self.insert1(key)
        print(f'Perform trial-segmentation of spike times for unit: {key["unit_id"]} and trial: {key["trial_id"]}')


@schema
class UnitPSTH(dj.Computed):
    definition = """  
    -> TrialSegmentedUnitSpikeTimes
    ---
    psth: longblob
    psth_time: longblob    
    """

    def make(self, key):
        return NotImplementedError



