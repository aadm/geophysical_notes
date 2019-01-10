"""
A Python module for reading/writing/manipulating
SEG-Y formatted files

segypy.read_segy                      : Read SEGY file
segypy.read_reel_header               : Get SEGY header
segypy.read_trace_header              : Get SEGY Trace header
segypy.read_all_trace_headers         : Get all SEGY Trace headers
segypy.get_default_segy_header        :
segypy.get_default_segy_trace_headers :
segypy.write_segy                     : Write a data to a SEGY file
segypy.write_segy_structure           : Writes a segpy data structure to a SEGY file
segypy.read_binary_value              : Get a value from a binary string
segypy.ibm2ieee                       : Convert IBM floats to IEEE
segypy.version                        : The version of SegyPY
segypy.verbose                        : Amount of verbose information to the screen
"""
#
# segpy : A Python module for reading and writing SEG-Y formatted data
#
# Forked by Robert Smallshire from the original segypy by
#
# (C) Thomas Mejer Hansen, 2005-2006
#
# with contributions from Pete Forman and Andrew Squelch 2007
#
# ++ modified by aadm, january 2015, to make it a single file
# ++ segypy.py now includes header_definition.py, ibm_float.py, revisions.py, trace_header_definition.py

import os

import sys
import struct
import logging

from numpy import (transpose, reshape, zeros, arange)

# from revisions import canonicalize_revision
SEGY_REVISION_0 = 0x0000
SEGY_REVISION_1 = 0x0100

VARIANTS = {SEGY_REVISION_0: SEGY_REVISION_0,
            SEGY_REVISION_1: SEGY_REVISION_1,
            100: SEGY_REVISION_1,
            256: SEGY_REVISION_1}

class SegYRevisionError(Exception):
    pass

def canonicalize_revision(revision):
    """Canonicalize a SEG Y revision.

    Various SEG Y revisions are seen in the wild; this function canonicalizes the supplies revision
    to either SEGY_REVISION_0 or SEGY_REVISION_1.

    Args:
        revision: Any object representing a SEG Y revision.

    Returns:
        Either SEGY_REVISION_0 or SEGY_REVISION_1.

    Raises:
        SegYRevisionError: If the revision is not known.
    """
    try:
        return VARIANTS[revision]
    except KeyError:
        raise SegYRevisionError("Unknown SEG Y Revision {!r}".format(revision))

# from header_definition import HEADER_DEF

HEADER_DEF = {"Job": {"pos": 3200, "type": "int32", "def": 0}}
HEADER_DEF["Line"] = {"pos": 3204, "type": "int32", "def": 0}
HEADER_DEF["Reel"] = {"pos": 3208, "type": "int32", "def": 0}
HEADER_DEF["DataTracePerEnsemble"] = {"pos": 3212, "type": "int16", "def": 0}
HEADER_DEF["AuxiliaryTracePerEnsemble"] = {"pos": 3214, "type": "int16", "def": 0}
HEADER_DEF["dt"] = {"pos": 3216, "type": "uint16", "def": 1000}
HEADER_DEF["dtOrig"] = {"pos": 3218, "type": "uint16", "def": 1000}
HEADER_DEF["ns"] = {"pos": 3220, "type": "uint16", "def": 0}
HEADER_DEF["nsOrig"] = {"pos": 3222, "type": "uint16", "def": 0}
HEADER_DEF["DataSampleFormat"] = {"pos": 3224, "type": "int16", "def": 5}
HEADER_DEF["DataSampleFormat"]["descr"] = {SEGY_REVISION_0: {
    1: "IBM Float",
    2: "32 bit Integer",
    3: "16 bit Integer",
    8: "8 bit Integer"}}

HEADER_DEF["DataSampleFormat"]["descr"][SEGY_REVISION_1] = {
    1: "IBM Float",
    2: "32 bit Integer",
    3: "16 bit Integer",
    5: "IEEE",
    8: "8 bit Integer"}

HEADER_DEF["DataSampleFormat"]["bps"] = {SEGY_REVISION_0: {
    1: 4,
    2: 4,
    3: 2,
    8: 1}}
HEADER_DEF["DataSampleFormat"]["bps"][SEGY_REVISION_1] = {
    1: 4,
    2: 4,
    3: 2,
    5: 4,
    8: 1}
HEADER_DEF["DataSampleFormat"]["datatype"] = {SEGY_REVISION_0: {
    1: 'ibm',
    2: 'l',
    3: 'h',
    8: 'B'}}
HEADER_DEF["DataSampleFormat"]["datatype"][SEGY_REVISION_1] = {
    1: 'ibm',
    2: 'l',
    3: 'h',
    #    5: 'float',
    5: 'f',
    8: 'B'}

HEADER_DEF["EnsembleFold"] = {"pos": 3226, "type": "int16", "def": 0}
HEADER_DEF["TraceSorting"] = {"pos": 3228, "type": "int16", "def": 0}
HEADER_DEF["VerticalSumCode"] = {"pos": 3230, "type": "int16", "def": 0}
HEADER_DEF["SweepFrequencyEnd"] = {"pos": 3234, "type": "int16", "def": 0}
HEADER_DEF["SweepLength"] = {"pos": 3236, "type": "int16", "def": 0}
HEADER_DEF["SweepType"] = {"pos": 3238, "type": "int16", "def": 0}
HEADER_DEF["SweepChannel"] = {"pos": 3240, "type": "int16", "def": 0}
HEADER_DEF["SweepTaperLengthStart"] = {"pos": 3242, "type": "int16", "def": 0}
HEADER_DEF["SweepTaperLengthEnd"] = {"pos": 3244, "type": "int16", "def": 0}
HEADER_DEF["TaperType"] = {"pos": 3246, "type": "int16", "def": 0}
HEADER_DEF["CorrelatedDataTraces"] = {"pos": 3248, "type": "int16", "def": 0}
HEADER_DEF["BinaryGain"] = {"pos": 3250, "type": "int16", "def": 0}
HEADER_DEF["AmplitudeRecoveryMethod"] = {"pos": 3252, "type": "int16", "def": 0}
HEADER_DEF["MeasurementSystem"] = {"pos": 3254, "type": "int16", "def": 0}
HEADER_DEF["ImpulseSignalPolarity"] = {"pos": 3256, "type": "int16", "def": 0}
HEADER_DEF["VibratoryPolarityCode"] = {"pos": 3258, "type": "int16", "def": 0}
HEADER_DEF["Unassigned1"] = {"pos": 3260, "type": "int16", "n": 120, "def": 0}
HEADER_DEF["SegyFormatRevisionNumber"] = {
    "pos": 3500, "type": "uint16", "def": 100}
HEADER_DEF["FixedLengthTraceFlag"] = {"pos": 3502, "type": "uint16", "def": 0}
HEADER_DEF["NumberOfExtTextualHeaders"] = {"pos": 3504, "type": "uint16", "def": 0}
HEADER_DEF["Unassigned2"] = {"pos": 3506, "type": "int16", "n": 47, "def": 0}

# from trace_header_definition import TRACE_HEADER_DEF

TRACE_HEADER_DEF = {"TraceSequenceLine": {"pos": 0, "type": "int32"}}
TRACE_HEADER_DEF["TraceSequenceFile"] = {"pos": 4, "type": "int32"}
TRACE_HEADER_DEF["FieldRecord"] = {"pos": 8, "type": "int32"}
TRACE_HEADER_DEF["TraceNumber"] = {"pos": 12, "type": "int32"}
TRACE_HEADER_DEF["EnergySourcePoint"] = {"pos": 16, "type": "int32"}
TRACE_HEADER_DEF["cdp"] = {"pos": 20, "type": "int32"}
TRACE_HEADER_DEF["cdpTrace"] = {"pos": 24, "type": "int32"}
TRACE_HEADER_DEF["TraceIdentificationCode"] = {"pos": 28, "type": "int16"}
TRACE_HEADER_DEF["TraceIdentificationCode"]["descr"] = {SEGY_REVISION_0: {
    1: "Seismic data",
    2: "Dead",
    3: "Dummy",
    4: "Time Break",
    5: "Uphole",
    6: "Sweep",
    7: "Timing",
    8: "Water Break"}}
TRACE_HEADER_DEF["TraceIdentificationCode"]["descr"][SEGY_REVISION_1] = {
    -1: "Other",
    0: "Unknown",
    1: "Seismic data",
    2: "Dead",
    3: "Dummy",
    4: "Time break",
    5: "Uphole",
    6: "Sweep",
    7: "Timing",
    8: "Waterbreak",
    9: "Near-field gun signature",
    10: "Far-field gun signature",
    11: "Seismic pressure sensor",
    12: "Multicomponent seismic sensor - Vertical component",
    13: "Multicomponent seismic sensor - Cross-line component",
    14: "Multicomponent seismic sensor - In-line component",
    15: "Rotated multicomponent seismic sensor - Vertical component",
    16: "Rotated multicomponent seismic sensor - Transverse component",
    17: "Rotated multicomponent seismic sensor - Radial component",
    18: "Vibrator reaction mass",
    19: "Vibrator baseplate",
    20: "Vibrator estimated ground force",
    21: "Vibrator reference",
    22: "Time-velocity pairs"}
TRACE_HEADER_DEF["NSummedTraces"] = {"pos": 30, "type": "int16"}
TRACE_HEADER_DEF["NStackedTraces"] = {"pos": 32, "type": "int16"}
TRACE_HEADER_DEF["DataUse"] = {"pos": 34, "type": "int16"}
TRACE_HEADER_DEF["DataUse"]["descr"] = {0: {
    1: "Production",
    2: "Test"}}
TRACE_HEADER_DEF["DataUse"]["descr"][1] = TRACE_HEADER_DEF["DataUse"]["descr"][0]
TRACE_HEADER_DEF["offset"] = {"pos": 36, "type": "int32"}
TRACE_HEADER_DEF["ReceiverGroupElevation"] = {"pos": 40, "type": "int32"}
TRACE_HEADER_DEF["SourceSurfaceElevation"] = {"pos": 44, "type": "int32"}
TRACE_HEADER_DEF["SourceDepth"] = {"pos": 48, "type": "int32"}
TRACE_HEADER_DEF["ReceiverDatumElevation"] = {"pos": 52, "type": "int32"}
TRACE_HEADER_DEF["SourceDatumElevation"] = {"pos": 56, "type": "int32"}
TRACE_HEADER_DEF["SourceWaterDepth"] = {"pos": 60, "type": "int32"}
TRACE_HEADER_DEF["GroupWaterDepth"] = {"pos": 64, "type": "int32"}
TRACE_HEADER_DEF["ElevationScalar"] = {"pos": 68, "type": "int16"}
TRACE_HEADER_DEF["SourceGroupScalar"] = {"pos": 70, "type": "int16"}
TRACE_HEADER_DEF["SourceX"] = {"pos": 72, "type": "int32"}
TRACE_HEADER_DEF["SourceY"] = {"pos": 76, "type": "int32"}
TRACE_HEADER_DEF["GroupX"] = {"pos": 80, "type": "int32"}
TRACE_HEADER_DEF["GroupY"] = {"pos": 84, "type": "int32"}
TRACE_HEADER_DEF["CoordinateUnits"] = {"pos": 88, "type": "int16"}
TRACE_HEADER_DEF["CoordinateUnits"]["descr"] = {SEGY_REVISION_0: {
    1: "Length (meters or feet)",
    2: "Seconds of arc"}}
TRACE_HEADER_DEF["CoordinateUnits"]["descr"][SEGY_REVISION_1] = {
    1: "Length (meters or feet)",
    2: "Seconds of arc",
    3: "Decimal degrees",
    4: "Degrees, minutes, seconds (DMS)"}
TRACE_HEADER_DEF["WeatheringVelocity"] = {"pos": 90, "type": "int16"}
TRACE_HEADER_DEF["SubWeatheringVelocity"] = {"pos": 92, "type": "int16"}
TRACE_HEADER_DEF["SourceUpholeTime"] = {"pos": 94, "type": "int16"}
TRACE_HEADER_DEF["GroupUpholeTime"] = {"pos": 96, "type": "int16"}
TRACE_HEADER_DEF["SourceStaticCorrection"] = {"pos": 98, "type": "int16"}
TRACE_HEADER_DEF["GroupStaticCorrection"] = {"pos": 100, "type": "int16"}
TRACE_HEADER_DEF["TotalStaticApplied"] = {"pos": 102, "type": "int16"}
TRACE_HEADER_DEF["LagTimeA"] = {"pos": 104, "type": "int16"}
TRACE_HEADER_DEF["LagTimeB"] = {"pos": 106, "type": "int16"}
TRACE_HEADER_DEF["DelayRecordingTime"] = {"pos": 108, "type": "int16"}
TRACE_HEADER_DEF["MuteTimeStart"] = {"pos": 110, "type": "int16"}
TRACE_HEADER_DEF["MuteTimeEND"] = {"pos": 112, "type": "int16"}
TRACE_HEADER_DEF["ns"] = {"pos": 114, "type": "uint16"}
TRACE_HEADER_DEF["dt"] = {"pos": 116, "type": "uint16"}
TRACE_HEADER_DEF["GainType"] = {"pos": 118, "type": "int16"}
TRACE_HEADER_DEF["GainType"]["descr"] = {SEGY_REVISION_0: {
    1: "Fixes",
    2: "Binary",
    3: "Floating point"}}
TRACE_HEADER_DEF["GainType"]["descr"][SEGY_REVISION_1] = TRACE_HEADER_DEF[
    "GainType"]["descr"][SEGY_REVISION_0]
TRACE_HEADER_DEF["InstrumentGainConstant"] = {"pos": 120, "type": "int16"}
TRACE_HEADER_DEF["InstrumentInitialGain"] = {"pos": 122, "type": "int16"}
TRACE_HEADER_DEF["Correlated"] = {"pos": 124, "type": "int16"}
TRACE_HEADER_DEF["Correlated"]["descr"] = {SEGY_REVISION_0: {
    1: "No",
    2: "Yes"}}
TRACE_HEADER_DEF["Correlated"]["descr"][SEGY_REVISION_1] = TRACE_HEADER_DEF[
    "Correlated"]["descr"][SEGY_REVISION_0]

TRACE_HEADER_DEF["SweepFrequencyStart"] = {"pos": 126, "type": "int16"}
TRACE_HEADER_DEF["SweepFrequencyEnd"] = {"pos": 128, "type": "int16"}
TRACE_HEADER_DEF["SweepLength"] = {"pos": 130, "type": "int16"}
TRACE_HEADER_DEF["SweepType"] = {"pos": 132, "type": "int16"}
TRACE_HEADER_DEF["SweepType"]["descr"] = {SEGY_REVISION_0: {
    1: "linear",
    2: "parabolic",
    3: "exponential",
    4: "other"}}
TRACE_HEADER_DEF["SweepType"]["descr"][SEGY_REVISION_1] = TRACE_HEADER_DEF[
    "SweepType"]["descr"][SEGY_REVISION_0]

TRACE_HEADER_DEF["SweepTraceTaperLengthStart"] = {"pos": 134, "type": "int16"}
TRACE_HEADER_DEF["SweepTraceTaperLengthEnd"] = {"pos": 136, "type": "int16"}
TRACE_HEADER_DEF["TaperType"] = {"pos": 138, "type": "int16"}
TRACE_HEADER_DEF["TaperType"]["descr"] = {SEGY_REVISION_0: {
    1: "linear",
    2: "cos2c",
    3: "other"}}
TRACE_HEADER_DEF["TaperType"]["descr"][SEGY_REVISION_1] = TRACE_HEADER_DEF[
    "TaperType"]["descr"][SEGY_REVISION_0]

TRACE_HEADER_DEF["AliasFilterFrequency"] = {"pos": 140, "type": "int16"}
TRACE_HEADER_DEF["AliasFilterSlope"] = {"pos": 142, "type": "int16"}
TRACE_HEADER_DEF["NotchFilterFrequency"] = {"pos": 144, "type": "int16"}
TRACE_HEADER_DEF["NotchFilterSlope"] = {"pos": 146, "type": "int16"}
TRACE_HEADER_DEF["LowCutFrequency"] = {"pos": 148, "type": "int16"}
TRACE_HEADER_DEF["HighCutFrequency"] = {"pos": 150, "type": "int16"}
TRACE_HEADER_DEF["LowCutSlope"] = {"pos": 152, "type": "int16"}
TRACE_HEADER_DEF["HighCutSlope"] = {"pos": 154, "type": "int16"}
TRACE_HEADER_DEF["YearDataRecorded"] = {"pos": 156, "type": "int16"}
TRACE_HEADER_DEF["DayOfYear"] = {"pos": 158, "type": "int16"}
TRACE_HEADER_DEF["HourOfDay"] = {"pos": 160, "type": "int16"}
TRACE_HEADER_DEF["MinuteOfHour"] = {"pos": 162, "type": "int16"}
TRACE_HEADER_DEF["SecondOfMinute"] = {"pos": 164, "type": "int16"}
TRACE_HEADER_DEF["TimeBaseCode"] = {"pos": 166, "type": "int16"}
TRACE_HEADER_DEF["TimeBaseCode"]["descr"] = {SEGY_REVISION_0: {
    1: "Local",
    2: "GMT",
    3: "Other"}}
TRACE_HEADER_DEF["TimeBaseCode"]["descr"][SEGY_REVISION_1] = {
    1: "Local",
    2: "GMT",
    3: "Other",
    4: "UTC"}
TRACE_HEADER_DEF["TraceWeightingFactor"] = {"pos": 168, "type": "int16"}
TRACE_HEADER_DEF["GeophoneGroupNumberRoll1"] = {"pos": 170, "type": "int16"}
TRACE_HEADER_DEF["GeophoneGroupNumberFirstTraceOrigField"] = {
    "pos": 172, "type": "int16"}
TRACE_HEADER_DEF["GeophoneGroupNumberLastTraceOrigField"] = {
    "pos": 174, "type": "int16"}
TRACE_HEADER_DEF["GapSize"] = {"pos": 176, "type": "int16"}
TRACE_HEADER_DEF["OverTravel"] = {"pos": 178, "type": "int16"}
TRACE_HEADER_DEF["OverTravel"]["descr"] = {SEGY_REVISION_0: {
    1: "down (or behind)",
    2: "up (or ahead)",
    3: "other"}}
TRACE_HEADER_DEF["OverTravel"]["descr"][SEGY_REVISION_1] = TRACE_HEADER_DEF[
    "OverTravel"]["descr"][SEGY_REVISION_0]


TRACE_HEADER_DEF["cdpX"] = {"pos": 180, "type": "int32"}
TRACE_HEADER_DEF["cdpY"] = {"pos": 184, "type": "int32"}
TRACE_HEADER_DEF["Inline3D"] = {"pos": 188, "type": "int32"}
TRACE_HEADER_DEF["Crossline3D"] = {"pos": 192, "type": "int32"}
TRACE_HEADER_DEF["ShotPoint"] = {"pos": 196, "type": "int32"}
TRACE_HEADER_DEF["ShotPointScalar"] = {"pos": 200, "type": "int16"}
TRACE_HEADER_DEF["TraceValueMeasurementUnit"] = {"pos": 202, "type": "int16"}
TRACE_HEADER_DEF["TraceValueMeasurementUnit"]["descr"] = {SEGY_REVISION_1: {
    -1: "Other",
    0: "Unknown (should be described in Data Sample "
       "Measurement Units Stanza)",
    1: "Pascal (Pa)",
    2: "Volts (V)",
    3: "Millivolts (v)",
    4: "Amperes (A)",
    5: "Meters (m)",
    6: "Meters Per Second (m/s)",
    7: "Meters Per Second squared (m/&s2)Other",
    8: "Newton (N)",
    9: "Watt (W)"}}
TRACE_HEADER_DEF["TransductionConstantMantissa"] = {"pos": 204, "type": "int32"}
TRACE_HEADER_DEF["TransductionConstantPower"] = {"pos": 208, "type": "int16"}
TRACE_HEADER_DEF["TransductionUnit"] = {"pos": 210, "type": "int16"}
TRACE_HEADER_DEF["TransductionUnit"]["descr"] = TRACE_HEADER_DEF[
    "TraceValueMeasurementUnit"]["descr"]
TRACE_HEADER_DEF["TraceIdentifier"] = {"pos": 212, "type": "int16"}
TRACE_HEADER_DEF["ScalarTraceHeader"] = {"pos": 214, "type": "int16"}
TRACE_HEADER_DEF["SourceType"] = {"pos": 216, "type": "int16"}
TRACE_HEADER_DEF["SourceType"]["descr"] = {SEGY_REVISION_1: {
    -1: "Other (should be described in Source Type/Orientation stanza)",
    0: "Unknown",
    1: "Vibratory - Vertical orientation",
    2: "Vibratory - Cross-line orientation",
    3: "Vibratory - In-line orientation",
    4: "Impulsive - Vertical orientation",
    5: "Impulsive - Cross-line orientation",
    6: "Impulsive - In-line orientation",
    7: "Distributed Impulsive - Vertical orientation",
    8: "Distributed Impulsive - Cross-line orientation",
    9: "Distributed Impulsive - In-line orientation"}}

TRACE_HEADER_DEF["SourceEnergyDirectionMantissa"] = {"pos": 218, "type": "int32"}
TRACE_HEADER_DEF["SourceEnergyDirectionExponent"] = {"pos": 222, "type": "int16"}
TRACE_HEADER_DEF["SourceMeasurementMantissa"] = {"pos": 224, "type": "int32"}
TRACE_HEADER_DEF["SourceMeasurementExponent"] = {"pos": 228, "type": "int16"}
TRACE_HEADER_DEF["SourceMeasurementUnit"] = {"pos": 230, "type": "int16"}
TRACE_HEADER_DEF["SourceMeasurementUnit"]["descr"] = {1: {
    -1: "Other (should be described in Source Measurement Unit stanza)",
    0: "Unknown",
    1: "Joule (J)",
    2: "Kilowatt (kW)",
    3: "Pascal (Pa)",
    4: "Bar (Bar)",
    5: "Newton (N)",
    6: "Kilograms (kg)"}}
TRACE_HEADER_DEF["UnassignedInt1"] = {"pos": 232, "type": "int32"}
TRACE_HEADER_DEF["UnassignedInt2"] = {"pos": 236, "type": "int32"}

# from ibm_float import ibm2ieee2
def ibm2ieee2(ibm_float):
    """
    ibm2ieee2(ibm_float)
    Used by permission
    (C) Secchi Angelo
    with thanks to Howard Lightstone and Anton Vredegoor.
    """
    dividend = float(16 ** 6)

    if ibm_float == 0:
        return 0.0
    istic, a, b, c = struct.unpack('>BBBB', ibm_float)
    if istic >= 128:
        sign = -1.0
        istic -= 128
    else:
        sign = 1.0
    mant = float(a << 16) + float(b << 8) + float(c)
    return sign * 16 ** (istic - 64) * (mant / dividend)

FORMAT = '%(asctime)-15s %(levelname)s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger('segpy.segypy')

version = '0.3.1'

REEL_HEADER_NUM_BYTES = 3600
TRACE_HEADER_NUM_BYTES = 240

l_char = struct.calcsize('c')
l_uchar = struct.calcsize('B')
l_float = struct.calcsize('f')


DATA_SAMPLE_FORMAT = {1: 'ibm',
                      2: 'l',
                      3: 'h',
                      5: 'f',
                      8: 'B'}

CTYPES = {'l': 'l', 'long':   'l', 'int32':  'l',
          'L': 'L', 'ulong':  'L', 'uint32': 'L',
          'h': 'h', 'short':  'h', 'int16':  'h',
          'H': 'H', 'ushort': 'H', 'uint16': 'H',
          'c': 'c', 'char':   'c',
          'B': 'B', 'uchar':  'B',
          'f': 'f', 'float':  'f',
          'ibm': 'ibm'}

# TODO This is redundant with data in the SH_def below
CTYPE_DESCRIPTION = {'ibm': 'IBM float',
                     'l':   '32 bit integer',
                     'h':   '16 bit integer',
                     'f':   'IEEE float',
                     'B':   '8 bit char'}


def size_in_bytes(ctype):
    if ctype == 'l' and struct.calcsize(ctype) == 8:
        return 4  # 64-bit issue?
    return struct.calcsize(ctype) if ctype != 'ibm' else struct.calcsize('f')


def get_default_segy_header(ntraces=100, ns=100):
    """
    header = getDefaultSegyHeader()
    """
    # TraceSequenceLine
    header = {'Job': {'pos': 3200, 'type': 'int32', 'def': 0}}

    for key in HEADER_DEF:
        header[key] = HEADER_DEF[key].get('def', 0)

    header['ntraces'] = ntraces
    header['ns'] = ns

    return header


def get_default_segy_trace_headers(ntraces=100, ns=100, dt=1000):
    """
    SH = getDefaultSegyTraceHeader()
    """
    # INITIALIZE DICTIONARY
    trace_header = {'TraceSequenceLine': {'pos': 0, 'type': 'int32'}}

    for key in TRACE_HEADER_DEF:
        trace_header[key] = zeros(ntraces)

    for a in range(ntraces):
        trace_header['TraceSequenceLine'][a] = a + 1
        trace_header['TraceSequenceFile'][a] = a + 1
        trace_header['FieldRecord'][a] = 1000
        trace_header['TraceNumber'][a] = a + 1
        trace_header['ns'][a] = ns
        trace_header['dt'][a] = dt
    return trace_header


def read_trace_header(f, reel_header, trace_header_name='cdp', endian='>'):
    """
    read_trace_header(reel_header, TraceHeaderName)
    """

    bps = get_byte_per_sample(reel_header)

    # MAKE SOME LOOKUP TABLE THAT HOLDS THE LOCATION OF HEADERS
    trace_header_pos = TRACE_HEADER_DEF[trace_header_name]['pos']

    # TODO: Be consistent between 'type' and 'format' here.
    trace_header_format = TRACE_HEADER_DEF[trace_header_name]['type']
    ntraces = reel_header['ntraces']
    trace_header_values = zeros(ntraces)
    binary_reader = create_binary_reader(f, trace_header_format, endian)
    start_pos = trace_header_pos + REEL_HEADER_NUM_BYTES
    stride = reel_header['ns'] * bps + TRACE_HEADER_NUM_BYTES
    end_pos = start_pos + (ntraces - 1) * stride + 1
    for i, pos in enumerate(xrange(start_pos, end_pos, stride)):
        trace_header_values[i] = binary_reader(pos)
    return trace_header_values


# TODO: Get the parameter ordering of reel_header and f to be consistent
def read_all_trace_headers(f, reel_header):
    trace_headers = {'filename': reel_header['filename']}

    logger.debug('read_all_trace_headers : '
                 'trying to get all segy trace headers')

    for key in TRACE_HEADER_DEF.keys():
        trace_header = read_trace_header(f, reel_header, key)
        trace_headers[key] = trace_header
        logger.info("read_all_trace_headers :  " + key)

    return trace_headers


def file_length(f):
    pos = f.tell()
    f.seek(0, os.SEEK_END)
    file_size = f.tell()
    f.seek(pos, os.SEEK_SET)
    return file_size


def _filename(f):
    return f.name if hasattr(f, 'name') else '<unknown>'


def read_segy(f, endian='>'):
    """
    data, header, trace_headers = read_reel_header(f)
    """

    file_size = file_length(f)

    logger.debug("read_segy : Length of data : {0}".format(file_size))

    reel_header = read_reel_header(f, endian)  # modified by A Squelch

    # GET TRACE
    index = REEL_HEADER_NUM_BYTES
    bytes_per_sample = get_byte_per_sample(reel_header)
    num_data = (file_size - REEL_HEADER_NUM_BYTES) / bytes_per_sample

    data, reel_header, trace_headers = read_traces(f,
                                                   reel_header,
                                                   num_data,
                                                   bytes_per_sample,
                                                   index,
                                                   endian)

    logger.debug("read_segy :  Read segy data")  # modified by A Squelch

    return data, reel_header, trace_headers


def read_traces(f,
                reel_header,
                num_data,
                bytes_per_sample,
                index,
                endian='>'):  # added by A Squelch
    """Read the trace data.

    values, SegyHeader, SegyTraceHeaders = read_traces(data,
                                                       reel_header,
                                                       num_data,
                                                       bytes_per_sample,
                                                       index)
    """

    # Calculate number of dummy samples needed to account for Trace Headers
    num_dummy_samples = TRACE_HEADER_NUM_BYTES / bytes_per_sample
    logger.debug("read_traces : num_dummy_samples = " + str(num_dummy_samples))

    # READ ALL SEGY TRACE HEADERS
    trace_headers = read_all_trace_headers(f, reel_header)

    logger.info("read_traces : Reading segy data")

    dsf = reel_header['DataSampleFormat']
    ctype = DATA_SAMPLE_FORMAT[dsf]
    description = CTYPE_DESCRIPTION[ctype]
    logger.debug("read_traces : Assuming DSF = {0}, {1}".format(
        dsf, description))
    values, _ = read_binary_value(f, index, ctype, endian, num_data)

    logger.debug("read_traces : - reshaping")
    values = reshape(values,
                     (reel_header['ntraces'],
                      reel_header['ns'] + num_dummy_samples))
    logger.debug("read_traces : - stripping header dummy data")
    values = values[:, num_dummy_samples:
                    (reel_header['ns'] + num_dummy_samples)]
    logger.debug("read_traces : - transposing")
    values = transpose(values)

    # SOMEONE NEEDS TO IMPLEMENT A NICER WAY DO DEAL WITH DSF = 8
    if reel_header['DataSampleFormat'] == 8:
        for i in arange(reel_header['ntraces']):
            for j in arange(reel_header['ns']):
                if values[i][j] > 128:
                    values[i][j] = values[i][j] - 256

    logger.debug("read_traces : Finished reading segy data")

    return values, reel_header, trace_headers


def read_reel_header(f, endian='>'):
    """
    reel_header = read_reel_header(file_handle)
    """
    filename = _filename(f)
    reel_header = {'filename': filename}
    for key in HEADER_DEF.keys():
        pos = HEADER_DEF[key]['pos']
        format = HEADER_DEF[key]['type']

        reel_header[key], index = read_binary_value(f, pos, format, endian)

        logger.debug(str(pos) + " " +
                     str(format) +
                     "  Reading " + key +
                     "=" + str(reel_header[key]))

    # SET NUMBER OF BYTES PER DATA SAMPLE
    bps = get_byte_per_sample(reel_header)

    file_size = file_length(f)
    ntraces = (file_size - REEL_HEADER_NUM_BYTES) / \
              (reel_header['ns'] * bps + TRACE_HEADER_NUM_BYTES)
    reel_header['ntraces'] = ntraces

    logger.debug('read_reel_header : successfully read ' + filename)

    return reel_header


def write_segy(filename, data, dt=1000, trace_header_in=None, header_in=None):
    """
    write_segy(filename, data, dt)

    Write SEGY

    See also read_segy

    (c) 2005, Thomas Mejer Hansen

    MAKE OPTIONAL INPUT FOR ALL SEGYHTRACEHEADER VALUES

    """
    if header_in is None:
        header_in = {}
    if trace_header_in is None:
        trace_header_in = {}

    logger.debug("write_segy : Trying to write " + filename)

    shape = data.shape
    ns = shape[0]
    ntraces = shape[1]
    print ntraces, ns

    header = get_default_segy_header(ntraces, ns)
    trace_header = get_default_segy_trace_headers(ntraces, ns, dt)

    # Add trace_header_in, if exists...
    for key in trace_header_in.keys():
        print key
        for a in range(ntraces):
            trace_header[key] = trace_header_in[key][a]

    # Add header_in, if exists...
    for key in header_in.keys():
        print key
        header[key] = header_in[key]

    write_segy_structure(filename, data, header, trace_header)


def write_segy_structure(filename,
                         data,
                         header,
                         trace_header,
                         endian='>'):  # modified by A Squelch
    """
    writeSegyStructure(filename, data, header, trace_header)

    Write SEGY file using SegyPy data structures

    See also readSegy

    (c) 2005, Thomas Mejer Hansen

    """

    logger.debug("writeSegyStructure : Trying to write " + filename)

    f = open(filename, 'wb')

    # VERBOSE INF
    revision = canonicalize_revision(header['SegyFormatRevisionNumber'])
    dsf = header['DataSampleFormat']

    try:  # block added by A Squelch
        data_descriptor = HEADER_DEF['DataSampleFormat']['descr'][revision][dsf]
    except KeyError:
        logging.critical("  An error has occurred interpreting a SEGY binary"
                         "header key")
        logging.critical("  Please check the Endian setting for this "
                         "file: {0}".format(header["filename"]))
        sys.exit()

    logger.debug("writeSegyStructure : SEG-Y revision = " + str(revision))
    logger.debug("writeSegyStructure : DataSampleFormat = " +
                 str(dsf) +
                 "(" + data_descriptor + ")")

    # WRITE SEGY HEADER

    for key in HEADER_DEF.keys():
        pos = HEADER_DEF[key]["pos"]
        format = HEADER_DEF[key]["type"]
        value = header[key]
        put_value(value, f, pos, format, endian)

    # SEGY TRACES
    ctype = HEADER_DEF['DataSampleFormat']['datatype'][revision][dsf]
    bps = HEADER_DEF['DataSampleFormat']['bps'][revision][dsf]

    sizeT = TRACE_HEADER_NUM_BYTES + header['ns'] * bps

    for itrace in range(header['ntraces']):
        index = REEL_HEADER_NUM_BYTES + itrace * sizeT
        logger.debug('Writing Trace #' +
                     str(itrace + 1) +
                     '/' + str(header['ntraces']))
        # WRITE SEGY TRACE HEADER
        for key in TRACE_HEADER_DEF.keys():
            pos = index + TRACE_HEADER_DEF[key]['pos']
            format = TRACE_HEADER_DEF[key]['type']
            value = trace_header[key][itrace]
            logger.debug(str(pos) + " " +
                         str(format) +
                         "  Writing " + key +
                         "=" + str(value))
            put_value(value, f, pos, format, endian)

        # Write Data
        cformat = endian + ctype
        for s in range(header['ns']):
            strVal = struct.pack(cformat, data[s, itrace])
            f.seek(index +
                   TRACE_HEADER_NUM_BYTES +
                   s * struct.calcsize(cformat))
            f.write(strVal)

    f.close()


def put_value(value, fileid, index, ctype='l', endian='>', number=1):
    """
    putValue(data, index, ctype, endian, number)
    """
    ctype = CTYPES[ctype]

    cformat = endian + ctype*number

    logger.debug('putValue : cformat :  ' + cformat + ' ctype = ' + ctype)

    str_val = struct.pack(cformat, value)
    fileid.seek(index)
    fileid.write(str_val)

    return 1


def create_binary_reader(f, ctype='l', endian='>'):
    """Create a unary callable which reads a given binary data type from a file.
    """
    ctype = CTYPES[ctype]
    size = size_in_bytes(ctype)

    cformat = endian + ctype

    def reader(index):
        f.seek(index, os.SEEK_SET)
        data = f.read(size)
        # TODO: Check the content of data before proceeding
        value = struct.unpack(cformat, data)
        return value[0]

    return reader


def read_binary_value(f, index, ctype='l', endian='>', number=1):
    """
    read_binary_value(data, index, ctype, endian, number)
    """

    ctype = CTYPES[ctype]

    size = size_in_bytes(ctype)

    cformat = endian + ctype * number

    logger.debug('read_binary_value : cformat :  ' + cformat)

    index_end = index + size * number

    f.seek(index, os.SEEK_SET)
    data = f.read(size * number)
    if ctype == 'ibm':
        # ASSUME IBM FLOAT DATA
        value = range(number)
        for i in arange(number):
            index_ibm = i * 4
            value[i] = ibm2ieee2(data[index_ibm: index_ibm + 4])
        # this returns an array as opposed to a tuple
    else:
        # TODO: Check the content of data before proceeding
        value = struct.unpack(cformat, data)

    if ctype == 'B':
        logger.warning('read_binary_value : '
                       'Inefficient use of 1 byte Integer...', 1)

    logger.debug('read_binary_value : ' +
                 'start = ' + str(index) +
                 ' size = ' + str(size) +
                 ' number = ' + str(number) +
                 ' value = ' + str(value) +
                 ' cformat = ' + str(cformat))

    if number == 1:
        return value[0], index_end
    else:
        return value, index_end


def get_byte_per_sample(header):
    revision = canonicalize_revision(header["SegyFormatRevisionNumber"])
    dsf = header["DataSampleFormat"]

    try:  # block added by A Squelch
        bps = HEADER_DEF["DataSampleFormat"]["bps"][revision][dsf]
    except KeyError:
        # TODO: This should not be a critical failure - should just convert
        # exception
        logging.critical("  An error has occurred interpreting a SEGY "
                         "binary header key")
        logging.critical("Please check the Endian setting for "
                         "this file: {0}".format(header["filename"]))
        sys.exit()

    logger.debug("getBytePerSample :  bps = " + str(bps))

    return bps


def main():
    filename = r'C:\Users\rjs\opendtectroot\Blake_Ridge_Hydrates_3D' \
               r'\stack_final_scaled50_int8.sgy'
    with open(filename, 'rb') as segy:
        data, header, trace_header = read_segy(segy)


if __name__ == '__main__':
    main()
