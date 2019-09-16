#!/usr/bin/env python3
import os

import sys
from datetime import datetime
from dateutil.tz import tzlocal
import pytz
import re
import numpy as np
import pandas as pd
import warnings
import tqdm

from pipeline import (reference, subject, acquisition, analysis,
                      extracellular, behavior, utilities)
import pynwb
from pynwb import NWBFile, NWBHDF5IO

warnings.filterwarnings('ignore', module='pynwb')

# ============================== SET CONSTANTS ==========================================
default_nwb_output_dir = os.path.join('data', 'NWB 2.0')
zero_zero_time = datetime.strptime('00:00:00', '%H:%M:%S').time()  # no precise time available
hardware_filter = 'Bandpass filtered 300-6K Hz'
institution = 'Janelia Research Campus'


def export_to_nwb(session_key, nwb_output_dir=default_nwb_output_dir, save=False, overwrite=True):
    this_session = (acquisition.Session & session_key).fetch1()

    identifier = '_'.join([this_session['subject_id'],
                           this_session['session_time'].strftime('%Y-%m-%d'),
                           str(this_session['session_id'])])
    # =============== General ====================
    # -- NWB file - a NWB2.0 file for each session
    nwbfile = NWBFile(
        session_description=this_session['session_note'],
        identifier=identifier,
        session_start_time=this_session['session_time'],
        file_create_date=datetime.now(tzlocal()),
        experimenter='; '.join((acquisition.Session.Experimenter & session_key).fetch('experimenter')),
        institution=institution,
        related_publications='https://doi.org/10.1038/s41586-018-0642-9')
    # -- subject
    subj = (subject.Subject & session_key).fetch1()
    nwbfile.subject = pynwb.file.Subject(
        subject_id=this_session['subject_id'],
        description=subj['subject_description'],
        genotype=' x '.join((subject.Subject.Allele & session_key).fetch('allele')),
        sex=subj['sex'],
        species=subj['species'])

    # =============== Extracellular ====================
    probe_insertion = ((extracellular.ProbeInsertion & session_key).fetch1()
                       if extracellular.ProbeInsertion & session_key
                       else None)
    if probe_insertion:
        probe = nwbfile.create_device(name = probe_insertion['probe_name'])
        electrode_group = nwbfile.create_electrode_group(
            name='; '.join([f'{probe_insertion["probe_name"]}: {str(probe_insertion["channel_counts"])}']),
            description = 'N/A',
            device = probe,
            location = '; '.join([f'{k}: {str(v)}' for k, v in
                                  (reference.BrainLocation & probe_insertion).fetch1().items()]))

        for chn in (reference.Probe.Channel & probe_insertion).fetch(as_dict=True):
            nwbfile.add_electrode(id=chn['channel_id'],
                                  group=electrode_group,
                                  filtering=hardware_filter,
                                  imp=-1.,
                                  x=0.0,  # not available from data
                                  y=0.0,  # not available from data
                                  z=0.0,  # not available from data
                                  location=electrode_group.location)

        # --- unit spike times ---
        nwbfile.add_unit_column(name='depth', description='depth this unit')
        nwbfile.add_unit_column(name='quality', description='spike width of this unit')
        nwbfile.add_unit_column(name='cell_type', description='cell type (e.g. wide width, narrow width spiking)')

        for unit in (extracellular.UnitSpikeTimes & probe_insertion).fetch(as_dict=True):
            # make an electrode table region (which electrode(s) is this unit coming from)
            nwbfile.add_unit(id=unit['unit_id'],
                             electrodes=(unit['channel_id']
                                         if isinstance(unit['channel_id'], np.ndarray) else [unit['channel_id']]),
                             depth=unit['unit_depth'],
                             quality=unit['unit_quality'],
                             cell_type=unit['unit_cell_type'],
                             spike_times=unit['spike_times'])

    # =============== Behavior ====================
    behavior_data = ((behavior.LickTimes & session_key).fetch1()
                     if behavior.LickTimes & session_key
                     else None)
    if behavior_data:
        behav_acq = pynwb.behavior.BehavioralEvents(name = 'lick_times')
        nwbfile.add_acquisition(behav_acq)
        [behavior_data.pop(k) for k in behavior.LickTimes.primary_key]
        for b_k, b_v in behavior_data.items():
            behav_acq.create_timeseries(name=b_k,
                                        unit='a.u.',
                                        conversion=1.0,
                                        data=np.full_like(b_v, 1),
                                        timestamps=b_v)

    # =============== TrialSet ====================
    # NWB 'trial' (of type dynamic table) by default comes with three mandatory attributes:
    #                                                                       'id', 'start_time' and 'stop_time'.
    # Other trial-related information needs to be added in to the trial-table as additional columns (with column name
    # and column description)
    if acquisition.TrialSet & session_key:
        # Get trial descriptors from TrialSet.Trial and TrialStimInfo
        trial_columns = [{'name': tag,
                          'description': re.sub('\s+:|\s+', ' ', re.search(
                              f'(?<={tag})(.*)', str(acquisition.TrialSet.Trial.heading)).group())}
                         for tag in acquisition.TrialSet.Trial.fetch(as_dict=True, limit=1)[0].keys()
                         if tag not in acquisition.TrialSet.Trial.primary_key + ['start_time', 'stop_time']]

        # Trial Events - discard 'trial_start' and 'trial_stop' as we already have start_time and stop_time
        trial_events = set(((acquisition.TrialSet.EventTime & session_key)
                            - [{'trial_event': 'trial_start'}, {'trial_event': 'trial_stop'}]).fetch('trial_event'))
        event_names = [{'name': e, 'description': d}
                       for e, d in zip(*(reference.ExperimentalEvent & [{'event': k}
                                                                        for k in trial_events]).fetch('event',
                                                                                                      'description'))]
        # Add new table columns to nwb trial-table for trial-label
        for c in trial_columns + event_names:
            nwbfile.add_trial_column(**c)

        # Add entry to the trial-table
        for trial in (acquisition.TrialSet.Trial & session_key).fetch(as_dict=True):
            events = dict(zip(*(acquisition.TrialSet.EventTime & trial
                                & [{'trial_event': e} for e in trial_events]).fetch('trial_event', 'event_time')))
            trial_tag_value = {**trial, **events, 'stop_time': np.nan}  # No stop_time available for this dataset

            trial_tag_value['id'] = trial_tag_value['trial_id']  # rename 'trial_id' to 'id'
            # convert None to np.nan since nwb fields does not take None
            for k, v in trial_tag_value.items():
                trial_tag_value[k] = v if v is not None else np.nan
            [trial_tag_value.pop(k) for k in acquisition.TrialSet.Trial.primary_key]
            nwbfile.add_trial(**trial_tag_value)

    # trial-aligned unit PSTH
    # create unit and trial "dynamic-table-region"
    unit_regions = {u: nwbfile.units.create_region(name='', region=[no], description='')
                    for no, u in enumerate(nwbfile.units.id.data)}
    trial_regions = {u: nwbfile.trials.create_region(name = '', region = [no], description = '')
                     for no, u in enumerate(nwbfile.trials.id.data)}

    PSTH = pynwb.core.DynamicTable(name='PSTH', description='trial-aligned unit PSTH')
    PSTH.add_column(name='unit_id', description='unit_id - link to the units table')
    PSTH.add_column(name='trial_id', description='trial_id - link to the trial table')
    PSTH.add_column(name='psth', description='trial-aligned unit PSTH')
    PSTH.add_column(name='psth_time', description = 'timestamps of the PSTH')

    for unit_id, trial_id, psth, psth_time in zip(*(extracellular.PSTH & session_key).fetch(
            'unit_id', 'trial_id', 'psth', 'psth_time')):
        PSTH.add_row({'unit_id': unit_regions[unit_id],
                      'trial_id': trial_regions[trial_id],
                      'psth': psth,
                      'psth_time': psth_time})

    psth_mod = nwbfile.create_processing_module(name='ecephys', description='trial-aligned unit PSTH')
    psth_mod.add_data_interface(PSTH)

    # =============== Write NWB 2.0 file ===============
    if save:
        save_file_name = ''.join([nwbfile.identifier, '.nwb'])
        if not os.path.exists(nwb_output_dir):
            os.makedirs(nwb_output_dir)
        if not overwrite and os.path.exists(os.path.join(nwb_output_dir, save_file_name)):
            return nwbfile
        with NWBHDF5IO(os.path.join(nwb_output_dir, save_file_name), mode = 'w') as io:
            io.write(nwbfile)
            print(f'Write NWB 2.0 file: {save_file_name}')

    return nwbfile


# ============================== EXPORT ALL ==========================================

if __name__ == '__main__':
    if len(sys.argv) > 1:
        nwb_outdir = sys.argv[1]
    else:
        nwb_outdir = default_nwb_output_dir

    for skey in acquisition.Session.fetch('KEY'):
        export_to_nwb(skey, nwb_output_dir=nwb_outdir, save=True)

