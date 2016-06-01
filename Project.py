#!/usr/bin/env python2
##################################################
# GNU Radio Python Flow Graph
# Title: FM Stereo Receiver
# Author: Brett Nicholas
# Generated: Tue May 31 22:32:12 2016
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from gnuradio import analog
from gnuradio import audio
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from gnuradio.wxgui import fftsink2
from gnuradio.wxgui import forms
from gnuradio.wxgui import scopesink2
from gnuradio.wxgui import waterfallsink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import wx


class Project(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="FM Stereo Receiver")

        ##################################################
        # Variables
        ##################################################
        self.smux_filt_samprate = smux_filt_samprate = 256e3
        self.smux_decim = smux_decim = 8
        self.samp_rate = samp_rate = 2.048e6
        self.right_gain = right_gain = 3
        self.left_gain = left_gain = 3
        self.bpf_base = bpf_base = 23e3
        self.RF_Gain = RF_Gain = 45
        self.CF = CF = 99.3e6

        ##################################################
        # Blocks
        ##################################################
        self._samp_rate_text_box = forms.text_box(
        	parent=self.GetWin(),
        	value=self.samp_rate,
        	callback=self.set_samp_rate,
        	label="Sample Rate: 1.024M, 1.4M, 1.8M, 1.92M, 2.048M, 2.4M & 2. 56M",
        	converter=forms.float_converter(),
        )
        self.GridAdd(self._samp_rate_text_box, 1, 0, 1, 1)
        _right_gain_sizer = wx.BoxSizer(wx.VERTICAL)
        self._right_gain_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_right_gain_sizer,
        	value=self.right_gain,
        	callback=self.set_right_gain,
        	label="R Audio Gain",
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._right_gain_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_right_gain_sizer,
        	value=self.right_gain,
        	callback=self.set_right_gain,
        	minimum=0,
        	maximum=5,
        	num_steps=100,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.GridAdd(_right_gain_sizer, 0, 1, 1, 1)
        self.notebook_0 = self.notebook_0 = wx.Notebook(self.GetWin(), style=wx.NB_TOP)
        self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "BB Spectrum")
        self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "Demod Spectrum")
        self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "Stereo Spectrum")
        self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "Stereo Signal")
        self.GridAdd(self.notebook_0, 2, 0, 1, 2)
        _left_gain_sizer = wx.BoxSizer(wx.VERTICAL)
        self._left_gain_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_left_gain_sizer,
        	value=self.left_gain,
        	callback=self.set_left_gain,
        	label="L Audio Gain",
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._left_gain_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_left_gain_sizer,
        	value=self.left_gain,
        	callback=self.set_left_gain,
        	minimum=0,
        	maximum=5,
        	num_steps=100,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.GridAdd(_left_gain_sizer, 0, 0, 1, 1)
        _RF_Gain_sizer = wx.BoxSizer(wx.VERTICAL)
        self._RF_Gain_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_RF_Gain_sizer,
        	value=self.RF_Gain,
        	callback=self.set_RF_Gain,
        	label="RF Gain",
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._RF_Gain_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_RF_Gain_sizer,
        	value=self.RF_Gain,
        	callback=self.set_RF_Gain,
        	minimum=0,
        	maximum=100,
        	num_steps=45,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.GridAdd(_RF_Gain_sizer, 1, 1, 1, 1)
        self.wxgui_waterfallsink2_0 = waterfallsink2.waterfall_sink_c(
        	self.notebook_0.GetPage(0).GetWin(),
        	baseband_freq=0,
        	dynamic_range=100,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=samp_rate,
        	fft_size=512,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="Baseband Waterfall",
        	size=(800,100),
        )
        self.notebook_0.GetPage(0).GridAdd(self.wxgui_waterfallsink2_0.win, 3, 0, 1, 2)
        self.wxgui_scopesink2_0 = scopesink2.scope_sink_f(
        	self.notebook_0.GetPage(3).GetWin(),
        	title="Scope Plot",
        	sample_rate=32e3,
        	v_scale=0,
        	v_offset=0,
        	t_scale=0,
        	ac_couple=False,
        	xy_mode=False,
        	num_inputs=2,
        	trig_mode=wxgui.TRIG_MODE_AUTO,
        	y_axis_label="Counts",
        	size=(800,500),
        )
        self.notebook_0.GetPage(3).Add(self.wxgui_scopesink2_0.win)
        self.wxgui_fftsink2_0_1 = fftsink2.fft_sink_f(
        	self.notebook_0.GetPage(2).GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=32e3,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="Difference FFT ",
        	peak_hold=False,
        )
        self.notebook_0.GetPage(2).Add(self.wxgui_fftsink2_0_1.win)
        self.wxgui_fftsink2_0_0_0 = fftsink2.fft_sink_f(
        	self.notebook_0.GetPage(1).GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=samp_rate/8,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="Demodulated FFT",
        	peak_hold=False,
        	size=(800,800),
        )
        self.notebook_0.GetPage(1).Add(self.wxgui_fftsink2_0_0_0.win)
        self.wxgui_fftsink2_0_0 = fftsink2.fft_sink_c(
        	self.notebook_0.GetPage(0).GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=samp_rate,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="Baseband FFT",
        	peak_hold=False,
        	size=(800,100),
        )
        self.notebook_0.GetPage(0).GridAdd(self.wxgui_fftsink2_0_0.win, 2, 0, 1, 2)
        self.wxgui_fftsink2_0 = fftsink2.fft_sink_f(
        	self.notebook_0.GetPage(2).GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=32e3,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="Sum FFT",
        	peak_hold=False,
        )
        self.notebook_0.GetPage(2).Add(self.wxgui_fftsink2_0.win)
        self.rfgain = blocks.multiply_const_vcc((RF_Gain, ))
        self.low_pass_filter_1_0 = filter.fir_filter_fff(smux_decim, firdes.low_pass(
        	1, smux_filt_samprate, 15e3, 500, firdes.WIN_HAMMING, 1))
        self.low_pass_filter_0 = filter.fir_filter_ccf(2, firdes.low_pass(
        	2, samp_rate/4, 100e3, 500, firdes.WIN_KAISER, 6.76))
        self.iir_filter_xxx_0 = filter.iir_filter_ccf((-0.00266, 0.00504, -0.00309, -0.00136, 0.00663, -0.01052, 0.01103, -0.00731, 0.00016, 0.00800, -0.01396, 0.01490, -0.00971, -0.00035, 0.01173, -0.01979, 0.02054, -0.01240, -0.00273, 0.01960, -0.03122, 0.03124, -0.01669, -0.01017, 0.04137, -0.06448, 0.06476, -0.02634, -0.07449, 0.33571, -0.00000, -0.33571, 0.07449, 0.02634, -0.06476, 0.06448, -0.04137, 0.01017, 0.01669, -0.03124, 0.03122, -0.01960, 0.00273, 0.01240, -0.02054, 0.01979, -0.01173, 0.00035, 0.00971, -0.01490, 0.01396, -0.00800, -0.00016, 0.00731, -0.01103, 0.01052, -0.00663, 0.00136, 0.00309, -0.00504, 0.00266
        ), (1 , ), False)
        self.fir_filter_xxx_0_0 = filter.fir_filter_ccf(4, (1,1,1,1))
        self.fir_filter_xxx_0_0.declare_sample_delay(0)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate,True)
        self.blocks_sub_xx_0 = blocks.sub_ff(1)
        self.blocks_multiply_xx_1_0 = blocks.multiply_vcc(1)
        self.blocks_multiply_xx_0 = blocks.multiply_vff(1)
        self.blocks_multiply_const_vxx_0_0 = blocks.multiply_const_vff((right_gain, ))
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_vff((left_gain, ))
        self.blocks_multiply_conjugate_cc_0 = blocks.multiply_conjugate_cc(1)
        self.blocks_float_to_complex_0 = blocks.float_to_complex(1)
        self.blocks_file_source_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, "/Users/bretttt/iCloud_drive/16S/engs110/project/radio_dat/IQ_Data_STEREO1", True)
        self.blocks_divide_xx_1 = blocks.divide_cc(1)
        self.blocks_delay_2 = blocks.delay(gr.sizeof_gr_complex*1, 30)
        self.blocks_complex_to_real_0 = blocks.complex_to_real(1)
        self.blocks_complex_to_mag_0 = blocks.complex_to_mag(1)
        self.blocks_complex_to_imag_0 = blocks.complex_to_imag(1)
        self.blocks_add_xx_0 = blocks.add_vff(1)
        self.blocks_add_const_vxx_0 = blocks.add_const_vcc((0.1, ))
        self.baseband_LPF = filter.fir_filter_fff(smux_decim, firdes.low_pass(
        	1, smux_filt_samprate, 15e3, 500, firdes.WIN_KAISER, 6.76))
        self.band_pass_filter_0_0_0 = filter.fir_filter_fcc(1, firdes.complex_band_pass(
        	1, smux_filt_samprate, 18000, 20000, 1000, firdes.WIN_KAISER, 1))
        self.band_pass_filter_0 = filter.fir_filter_fff(1, firdes.band_pass(
        	1, smux_filt_samprate, bpf_base, bpf_base+30e3, 500, firdes.WIN_KAISER, 6.76))
        self.audio_sink_0_0_0_0 = audio.sink(32000, "", True)
        self.analog_pll_refout_cc_0_0 = analog.pll_refout_cc(3.14/100, 0.152*3.14, 0.144*3.14)
        self.analog_fm_deemph_0_0 = analog.fm_deemph(fs=samp_rate/8, tau=75e-6)
        self.analog_fm_deemph_0 = analog.fm_deemph(fs=samp_rate/8, tau=75e-6)
        self.analog_const_source_x_0 = analog.sig_source_f(0, analog.GR_CONST_WAVE, 0, 0, 0)
        _CF_sizer = wx.BoxSizer(wx.VERTICAL)
        self._CF_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_CF_sizer,
        	value=self.CF,
        	callback=self.set_CF,
        	label="Center Frequency",
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._CF_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_CF_sizer,
        	value=self.CF,
        	callback=self.set_CF,
        	minimum=80e6,
        	maximum=108e6,
        	num_steps=280,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.GridAdd(_CF_sizer, 3, 0, 1, 2)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_const_source_x_0, 0), (self.blocks_float_to_complex_0, 1))    
        self.connect((self.analog_fm_deemph_0, 0), (self.audio_sink_0_0_0_0, 0))    
        self.connect((self.analog_fm_deemph_0_0, 0), (self.audio_sink_0_0_0_0, 1))    
        self.connect((self.analog_pll_refout_cc_0_0, 0), (self.blocks_multiply_xx_1_0, 0))    
        self.connect((self.analog_pll_refout_cc_0_0, 0), (self.blocks_multiply_xx_1_0, 1))    
        self.connect((self.band_pass_filter_0, 0), (self.blocks_multiply_xx_0, 0))    
        self.connect((self.band_pass_filter_0_0_0, 0), (self.analog_pll_refout_cc_0_0, 0))    
        self.connect((self.baseband_LPF, 0), (self.blocks_add_xx_0, 0))    
        self.connect((self.baseband_LPF, 0), (self.blocks_sub_xx_0, 0))    
        self.connect((self.blocks_add_const_vxx_0, 0), (self.blocks_divide_xx_1, 1))    
        self.connect((self.blocks_add_xx_0, 0), (self.blocks_multiply_const_vxx_0, 0))    
        self.connect((self.blocks_complex_to_imag_0, 0), (self.band_pass_filter_0, 0))    
        self.connect((self.blocks_complex_to_imag_0, 0), (self.band_pass_filter_0_0_0, 0))    
        self.connect((self.blocks_complex_to_imag_0, 0), (self.baseband_LPF, 0))    
        self.connect((self.blocks_complex_to_imag_0, 0), (self.wxgui_fftsink2_0_0_0, 0))    
        self.connect((self.blocks_complex_to_mag_0, 0), (self.blocks_float_to_complex_0, 0))    
        self.connect((self.blocks_complex_to_real_0, 0), (self.blocks_multiply_xx_0, 1))    
        self.connect((self.blocks_delay_2, 0), (self.blocks_multiply_conjugate_cc_0, 1))    
        self.connect((self.blocks_divide_xx_1, 0), (self.blocks_delay_2, 0))    
        self.connect((self.blocks_divide_xx_1, 0), (self.iir_filter_xxx_0, 0))    
        self.connect((self.blocks_file_source_0_0, 0), (self.rfgain, 0))    
        self.connect((self.blocks_float_to_complex_0, 0), (self.blocks_add_const_vxx_0, 0))    
        self.connect((self.blocks_multiply_conjugate_cc_0, 0), (self.blocks_complex_to_imag_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.analog_fm_deemph_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.wxgui_fftsink2_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.wxgui_scopesink2_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0_0, 0), (self.analog_fm_deemph_0_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0_0, 0), (self.wxgui_fftsink2_0_1, 0))    
        self.connect((self.blocks_multiply_const_vxx_0_0, 0), (self.wxgui_scopesink2_0, 1))    
        self.connect((self.blocks_multiply_xx_0, 0), (self.low_pass_filter_1_0, 0))    
        self.connect((self.blocks_multiply_xx_1_0, 0), (self.blocks_complex_to_real_0, 0))    
        self.connect((self.blocks_sub_xx_0, 0), (self.blocks_multiply_const_vxx_0_0, 0))    
        self.connect((self.blocks_throttle_0, 0), (self.fir_filter_xxx_0_0, 0))    
        self.connect((self.blocks_throttle_0, 0), (self.wxgui_fftsink2_0_0, 0))    
        self.connect((self.blocks_throttle_0, 0), (self.wxgui_waterfallsink2_0, 0))    
        self.connect((self.fir_filter_xxx_0_0, 0), (self.low_pass_filter_0, 0))    
        self.connect((self.iir_filter_xxx_0, 0), (self.blocks_multiply_conjugate_cc_0, 0))    
        self.connect((self.low_pass_filter_0, 0), (self.blocks_complex_to_mag_0, 0))    
        self.connect((self.low_pass_filter_0, 0), (self.blocks_divide_xx_1, 0))    
        self.connect((self.low_pass_filter_1_0, 0), (self.blocks_add_xx_0, 1))    
        self.connect((self.low_pass_filter_1_0, 0), (self.blocks_sub_xx_0, 1))    
        self.connect((self.rfgain, 0), (self.blocks_throttle_0, 0))    


    def get_smux_filt_samprate(self):
        return self.smux_filt_samprate

    def set_smux_filt_samprate(self, smux_filt_samprate):
        self.smux_filt_samprate = smux_filt_samprate
        self.band_pass_filter_0.set_taps(firdes.band_pass(1, self.smux_filt_samprate, self.bpf_base, self.bpf_base+30e3, 500, firdes.WIN_KAISER, 6.76))
        self.band_pass_filter_0_0_0.set_taps(firdes.complex_band_pass(1, self.smux_filt_samprate, 18000, 20000, 1000, firdes.WIN_KAISER, 1))
        self.baseband_LPF.set_taps(firdes.low_pass(1, self.smux_filt_samprate, 15e3, 500, firdes.WIN_KAISER, 6.76))
        self.low_pass_filter_1_0.set_taps(firdes.low_pass(1, self.smux_filt_samprate, 15e3, 500, firdes.WIN_HAMMING, 1))

    def get_smux_decim(self):
        return self.smux_decim

    def set_smux_decim(self, smux_decim):
        self.smux_decim = smux_decim

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self._samp_rate_text_box.set_value(self.samp_rate)
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)
        self.low_pass_filter_0.set_taps(firdes.low_pass(2, self.samp_rate/4, 100e3, 500, firdes.WIN_KAISER, 6.76))
        self.wxgui_fftsink2_0_0.set_sample_rate(self.samp_rate)
        self.wxgui_fftsink2_0_0_0.set_sample_rate(self.samp_rate/8)
        self.wxgui_waterfallsink2_0.set_sample_rate(self.samp_rate)

    def get_right_gain(self):
        return self.right_gain

    def set_right_gain(self, right_gain):
        self.right_gain = right_gain
        self._right_gain_slider.set_value(self.right_gain)
        self._right_gain_text_box.set_value(self.right_gain)
        self.blocks_multiply_const_vxx_0_0.set_k((self.right_gain, ))

    def get_left_gain(self):
        return self.left_gain

    def set_left_gain(self, left_gain):
        self.left_gain = left_gain
        self._left_gain_slider.set_value(self.left_gain)
        self._left_gain_text_box.set_value(self.left_gain)
        self.blocks_multiply_const_vxx_0.set_k((self.left_gain, ))

    def get_bpf_base(self):
        return self.bpf_base

    def set_bpf_base(self, bpf_base):
        self.bpf_base = bpf_base
        self.band_pass_filter_0.set_taps(firdes.band_pass(1, self.smux_filt_samprate, self.bpf_base, self.bpf_base+30e3, 500, firdes.WIN_KAISER, 6.76))

    def get_RF_Gain(self):
        return self.RF_Gain

    def set_RF_Gain(self, RF_Gain):
        self.RF_Gain = RF_Gain
        self._RF_Gain_slider.set_value(self.RF_Gain)
        self._RF_Gain_text_box.set_value(self.RF_Gain)
        self.rfgain.set_k((self.RF_Gain, ))

    def get_CF(self):
        return self.CF

    def set_CF(self, CF):
        self.CF = CF
        self._CF_slider.set_value(self.CF)
        self._CF_text_box.set_value(self.CF)


if __name__ == '__main__':
    parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
    (options, args) = parser.parse_args()
    tb = Project()
    tb.Start(True)
    tb.Wait()
