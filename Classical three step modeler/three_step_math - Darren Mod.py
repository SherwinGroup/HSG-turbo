# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 09:53:27 2016

@author: dreadnought

Changed it to add a graph showing the trace of the THz field with arrows
indicating the creation/recollision times.
"""

import sys
import json
from PyQt5 import QtWidgets
import pyqtgraph as pg
from three_step_first_ui import Ui_Form
import numpy as np
import scipy.interpolate as spi
import scipy.optimize as spo

import matplotlib.pyplot as plt

import InstsAndQt

q_e = 1.602e-19 # Charge of electron in Coulombs
m_e = 9.109e-31 # Mass of bare electron in kg
hbar = 1.055e-34 # Reduced Planck's constant in J*s

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
class Window(QtWidgets.QWidget):
    def __init__(self):
        """
        I have no real understanding of the things that are happening, but I know
        they're the best way.  
        """
        super(Window, self).__init__()
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.show()
        # self.data = Trajectories(0.063, 0.51, 0.54, 10, 10, 60)
        self.ui.graphicsView_pos_graph.plotItem.setLabel('bottom', text='Time', units='s')
        self.ui.graphicsView_pos_graph.plotItem.setLabel('left', text='Position', units='m')
        self.ui.graphicsView_pos_graph.plotItem.addLegend()
        self.e_pos_line = self.ui.graphicsView_pos_graph.plot(pen='r', name='Electron trajectory')
        self.h_pos_line = self.ui.graphicsView_pos_graph.plot(pen='b', name='Hole trajectory')


        # self.data = Trajectories(0.063, 0.51, 0.54, 10, 10, 60)
        self.ui.graphicsView_field_graph.plotItem.setLabel('bottom', text='Time', units='s')
        self.ui.graphicsView_field_graph.plotItem.setLabel('left', text='Electric Field', units='V/cm')
        self.e_field_line = self.ui.graphicsView_field_graph.plot(pen='k')
        self.ionization_time_arrow = pg.ArrowItem(angle=90, tipAngle=60, headLen=10, tailLen=25, tailWidth=2, brush='r', pen=None)
        self.recollision_time_arrow = pg.ArrowItem(angle=90, tipAngle=60, headLen=10, tailLen=25, tailWidth=2, brush='b', pen=None)
        self.ui.graphicsView_field_graph.plotItem.addItem(self.ionization_time_arrow)
        self.ui.graphicsView_field_graph.plotItem.addItem(self.recollision_time_arrow)


        # self.pos_line.setData(self.data.time_vals_ps, self.data.e_position_nm)
        self.ui.graphicsView_energy_graph.plotItem.setLabel('bottom', text='Time', units='s')
        self.ui.graphicsView_energy_graph.plotItem.setLabel('left', text='Energy', units='eV')
        self.ui.graphicsView_energy_graph.plotItem.addLegend()
        self.tot_energy_line = self.ui.graphicsView_energy_graph.plot(pen='k', name='Total kinetic energy')
        self.e_energy_line = self.ui.graphicsView_energy_graph.plot(pen='r', name='Electron kinetic energy')
        self.h_energy_line = self.ui.graphicsView_energy_graph.plot(pen='b', name='Hole kinetic energy')
        
        self.get_changes()
        # self.energy_line.setData(self.data.time_vals_ps, self.data.e_kenergy_meV)
        self.ui.lineEdit_m_e_eff.editingFinished.connect(self.get_changes)
        self.ui.lineEdit_m_h_eff.editingFinished.connect(self.get_changes)
        self.ui.lineEdit_thz_freq.editingFinished.connect(self.get_changes)
        self.ui.lineEdit_thz_field.editingFinished.connect(self.get_changes)
        self.ui.lineEdit_e_bind.editingFinished.connect(self.get_changes)
        self.ui.lineEdit_phase.editingFinished.connect(self.get_changes)

        self.ui.lineEdit_tot_ke_rec.editingFinished.connect(self.calculatePhaseFromEnergy)
        self.ui.lineEdit_approx_sb.editingFinished.connect(self.calculatePhaseFromSBNum)
        self.ui.lineEdit_ionization_time.editingFinished.connect(self.calculatePhaseFromIonizationTime)
        self.ui.lineEdit_recollision_time.editingFinished.connect(self.calculatePhaseFromRecollisionTime)

        self.ui.lineEdit_slew_rate_time_fs.editingFinished.connect(self.get_slew_rate)
        self.ui.lineEdit_splitting_meV.editingFinished.connect(self.get_slew_rate)

        self.ui.pushButton_chooseFolder.clicked.connect(self.choose_folder) # Choose the path
        self.ui.pushButton_save.clicked.connect(self.save_file) # add the file name to the path
        self.ui.lineEdit_file_path.editingFinished.connect(self.set_file_name) # Save it all!

    def get_changes(self):
        """
        This is the method that is used as the slot for all of those lineEdit
        boxes.  I'm not super clear on why they are they way they are, but it's
        working so far.
        """
        m_e_eff = float(self.ui.lineEdit_m_e_eff.text())
        m_h_eff = float(self.ui.lineEdit_m_h_eff.text())
        thz_freq = float(self.ui.lineEdit_thz_freq.text())
        thz_field = float(self.ui.lineEdit_thz_field.text())
        e_bind = float(self.ui.lineEdit_e_bind.text())
        phi = float(self.ui.lineEdit_phase.text())

        self.data = Trajectories(m_e_eff, m_h_eff, thz_freq, thz_field, e_bind, phi)

        self.e_pos_line.setData(self.data.time_vals_s, self.data.e_position_m)
        self.h_pos_line.setData(self.data.time_vals_s, self.data.h_position_m)

        self.tot_energy_line.setData(self.data.time_vals_s, self.data.total_kenergy_eV)
        self.e_energy_line.setData(self.data.time_vals_s, self.data.e_kenergy_eV)
        self.h_energy_line.setData(self.data.time_vals_s, self.data.h_kenergy_eV)

        self.data.make_parameters()
        self.ui.lineEdit_approx_sb.setText('{:.1f}'.format(self.data.approx_sb))
        self.ui.lineEdit_recollision_duration.setText('{:.3f}'.format(self.data.recollision_time_ps))
        self.ui.lineEdit_e_ke_rec.setText('{:.3f}'.format(self.data.e_kenergy_rec_meV))
        self.ui.lineEdit_h_ke_rec.setText('{:.3f}'.format(self.data.h_kenergy_rec_meV))
        self.ui.lineEdit_tot_ke_rec.setText('{:.3f}'.format(self.data.total_kenergy_rec_meV))
        self.ui.lineEdit_keldysh.setText('{:.3f}'.format(self.data.keldysh))
        self.ui.lineEdit_ponderomotive.setText('{:.3f}'.format(self.data.ponderomotive_meV))
        self.ui.lineEdit_max_sbs.setText('{:.3f}'.format(self.data.max_sb_order))

        time_vals = np.linspace(-1, 1, 200)/(self.data.thz_freq*2)
        sin_vals = -self.data.thz_field*np.sin(2*np.pi*self.data.thz_freq*time_vals)/100. # v/cm
        self.e_field_line.setData(time_vals, sin_vals)

        # ion_time = (self.data.phi-np.pi/2)/(2*np.pi*self.data.thz_freq)
        ion_time = self.data.ion_time
        ion_field = -self.data.thz_field*np.sin(2*np.pi*self.data.thz_freq*ion_time)/100. # v/cm
        self.ionization_time_arrow.setPos(ion_time, ion_field)
        # recol_time = ion_time+self.data.recollision_time_ps*1e-12
        recol_time = self.data.recol_time
        recol_field = -self.data.thz_field*np.sin(2*np.pi*self.data.thz_freq*recol_time)/100. # v/cm
        self.recollision_time_arrow.setPos(recol_time, recol_field)

        self.ui.lineEdit_ionization_time.setText('{:.3f}'.format(ion_time*1e15))
        self.ui.lineEdit_recollision_time.setText('{:.3f}'.format(recol_time * 1e15))

    def get_slew_rate(self):
        time = float(self.ui.lineEdit_slew_rate_time_fs.text())
        splitting = float(self.ui.lineEdit_splitting_meV.text())
        self.data.make_energy_slew_rate(time, splitting)
        self.ui.lineEdit_slew_rate.setText('{:.3f}'.format(self.data.energy_slew_rate_meV_fs))
        self.ui.lineEdit_diabatic_prob.setText('{:.3f}'.format(self.data.diabatic_probability))

    def calculatePhaseFromIonizationTime(self):
        thz_freq = float(self.ui.lineEdit_thz_freq.text())
        ion_time = float(self.ui.lineEdit_ionization_time.text())*1e-3

        phi = (2*np.pi*thz_freq * ion_time + np.pi/2)*180./np.pi

        self.ui.lineEdit_ionization_time.blockSignals(True)
        self.ui.lineEdit_phase.setText("{:.3f}".format(phi))
        self.ui.lineEdit_ionization_time.blockSignals(False)
        self.get_changes()

    def calculatePhaseFromRecollisionTime(self):
        recoltime = float(self.ui.lineEdit_recollision_time.text())*1e-3

        m_e_eff = float(self.ui.lineEdit_m_e_eff.text())
        m_h_eff = float(self.ui.lineEdit_m_h_eff.text())
        thz_freq = float(self.ui.lineEdit_thz_freq.text())
        thz_field = float(self.ui.lineEdit_thz_field.text())
        e_bind = float(self.ui.lineEdit_e_bind.text())

        def minir(phi):
            try:
                data = Trajectories(m_e_eff, m_h_eff, thz_freq, thz_field, e_bind, phi)
                ion_time = (phi*np.pi/180. - np.pi / 2) / (2 * np.pi * thz_freq)
            except IndexError:
                return np.inf
            return float(data.recollision_time_ps) + ion_time

        f = lambda ph: np.abs(recoltime-minir(ph))
        p = spo.minimize(f, 60, method="Nelder-Mead")

        self.ui.lineEdit_recollision_time.blockSignals(True)
        self.ui.lineEdit_phase.setText("{:.3f}".format(p.x[0]))
        self.ui.lineEdit_recollision_time.blockSignals(False)
        self.get_changes()

    def calculatePhaseFromEnergy(self):
        energy = float(self.ui.lineEdit_tot_ke_rec.text())

        m_e_eff = float(self.ui.lineEdit_m_e_eff.text())
        m_h_eff = float(self.ui.lineEdit_m_h_eff.text())
        thz_freq = float(self.ui.lineEdit_thz_freq.text())
        thz_field = float(self.ui.lineEdit_thz_field.text())
        e_bind = float(self.ui.lineEdit_e_bind.text())

        def minir(phi):
            try:
                data = Trajectories(m_e_eff, m_h_eff, thz_freq, thz_field, e_bind, phi)
            except IndexError:
                return np.inf
            return float(data.total_kenergy_rec_meV)

        f = lambda ph: np.abs(energy-minir(ph))

        p = spo.minimize(f, 60, method="Nelder-Mead")

        self.ui.lineEdit_tot_ke_rec.blockSignals(True)
        self.ui.lineEdit_phase.setText("{:.3f}".format(p.x[0]))
        self.ui.lineEdit_tot_ke_rec.blockSignals(False)
        self.get_changes()

    def calculatePhaseFromSBNum(self):
        sb = float(self.ui.lineEdit_approx_sb.text())

        m_e_eff = float(self.ui.lineEdit_m_e_eff.text())
        m_h_eff = float(self.ui.lineEdit_m_h_eff.text())
        thz_freq = float(self.ui.lineEdit_thz_freq.text())
        thz_field = float(self.ui.lineEdit_thz_field.text())
        e_bind = float(self.ui.lineEdit_e_bind.text())

        def minir(phi):
            try:
                data = Trajectories(m_e_eff, m_h_eff, thz_freq, thz_field, e_bind, phi)
            except IndexError:
                return np.inf
            return float(data.approx_sb)

        f = lambda ph: np.abs(sb-minir(ph))

        p = spo.minimize(f, 60, method="Nelder-Mead")

        # self.ui.lineEdit_tot_ke_rec.blockSignals(True)
        self.ui.lineEdit_phase.setText("{:.3f}".format(p.x[0]))
        # self.ui.lineEdit_tot_ke_rec.blockSignals(False)
        self.get_changes()

    def choose_folder(self):
        hint = 'Choose save directory'
        self.file_dir = str(QtGui.QFileDialog.getExistingDirectory(self, hint))
        self.ui.lineEdit_file_path.setText(self.file_dir)

    def set_file_name(self):
        self.file_name = str(self.ui.lineEdit_file_path.text())

    def save_file(self):
        header_string = '#Initial settings:\n'
        header_string += json.dumps(self.data.parts_dict, sort_keys=True, indent=4, separators=(',', ': '))
        header_string += '\nResults:\n'
        header_string += json.dumps(self.data.results_dict, sort_keys=True, indent=4, separators=(',', ': '))
        header_string = header_string.replace('\n', '\n#')
        header_string += '\nTime,Electric field,Position,Position,Kinetic energy,' \
                         'Kinetic energy,Kinetic energy,Full Field Time,Full Field'
        header_string += '\nps,kV/cm,nm,nm,meV,meV,meV,ps,kV/cm'
        header_string += '\n,,Electron,Hole,Electron,Hole,Total,,{}/{}'.format(
            self.data.thz_freq*1e-9, self.data.thz_field*1e-5
        )

        time_vals_full_field = np.linspace(-1, 1, len(self.data.time_vals_s)) / (self.data.thz_freq * 2)
        sin_vals_full_field = -self.data.thz_field * np.sin(
            2 * np.pi * self.data.thz_freq * time_vals_full_field) / 100.  # v/cm


        save_stuff = np.column_stack((self.data.time_vals_s * 1e12, self.data.driving_field * self.data.thz_field / 1e5, self.data.e_position_m * 1e9,
                                      self.data.h_position_m * 1e9, self.data.e_kenergy_eV * 1e3, self.data.h_kenergy_eV * 1e3,
                                      self.data.total_kenergy_eV * 1e3,
                                      time_vals_full_field*1e12,
                                      sin_vals_full_field*1e-3))
        np.savetxt(self.file_name, save_stuff, comments='', header=header_string, delimiter=',')

class Trajectories(object):
    """
    This guy includes all the information about the three step Trajectories
    """
    def __init__(self, m_e_eff, m_h_eff, thz_freq, thz_field, e_bind, phi):
        """
        I have no idea what to include in here.
        """

        phi = float(phi)
        self.m_e_eff = m_e_eff * m_e # in kg
        self.m_h_eff = m_h_eff * m_e # in kg
        self.reduced_mass = 1 / (1 / self.m_e_eff + 1 / self.m_h_eff) # in kg
        self.thz_freq = thz_freq * 1e12 # in Hz
        self.thz_field = thz_field * 1e5 # in V/m
        self.e_bind = e_bind * q_e * 1e-3 # in J
        self.phi = phi * np.pi / 180 # in rad

        self.parts_dict = {'eff. electron mass': '{:.4f}'.format(m_e_eff),
                           'eff. hole mass': '{:.4f}'.format(m_h_eff),
                           'eff. reduced mass': '{:.4f}'.format(self.reduced_mass / m_e),
                           'thz field (kV/cm)': '{:.2f}'.format(thz_field),
                           'thz frequency (THz)': '{:.3f}'.format(thz_freq),
                           'exciton binding energy (eV)': '{:.3f}'.format(e_bind),
                           'tunneling phase (deg)': '{:.3f}'.format(phi)}
        self.results_dict = {}
        self.update_stuff()

    def update_stuff(self):
        """
        Updates/creates all of the relevant plots and parameters
        """
        self.make_parameters()

        period = 1 / self.thz_freq
        time_s = np.linspace(0, 1 * period, num=10000) # in seconds
        driving_field = np.cos(2 * np.pi * self.thz_freq * time_s + self.phi) # this is dimless
        position = -driving_field - 2 * np.pi * self.thz_freq * time_s * np.sin(self.phi) + np.cos(self.phi)

        recollision_index = np.where(position < 0)[0][0]
        self.recollision_time_ps = time_s[recollision_index] * 1e12

        self.time_vals_s = np.linspace(0, time_s[recollision_index], num=1000)
        self.driving_field = np.cos(2 * np.pi * self.thz_freq * self.time_vals_s + self.phi) # this is dimless

        position = -self.driving_field - 2 * np.pi * self.thz_freq * self.time_vals_s * np.sin(self.phi) + np.cos(self.phi)
        e_pos_scale_m = q_e * self.thz_field / (self.m_e_eff * (2 * np.pi * self.thz_freq)**2)
        h_pos_scale_m = -q_e * self.thz_field / (self.m_h_eff * (2 * np.pi * self.thz_freq)**2) # Yes, I made electron go positive
        self.e_position_m = e_pos_scale_m * position
        self.h_position_m = h_pos_scale_m * position
        
        e_kenergy_scale = 1 / (2 * self.m_e_eff) * (q_e * self.thz_field / (2 * np.pi * self.thz_freq))**2 / q_e
        h_kenergy_scale = 1 / (2 * self.m_h_eff) * (q_e * self.thz_field / (2 * np.pi * self.thz_freq))**2 / q_e
        energy = (np.sin(2 * np.pi * self.thz_freq * self.time_vals_s + self.phi) - np.sin(self.phi))**2
        self.e_kenergy_eV = e_kenergy_scale * energy
        self.h_kenergy_eV = h_kenergy_scale * energy
        self.total_kenergy_eV = self.e_kenergy_eV + self.h_kenergy_eV
        self.e_kenergy_rec_meV = self.e_kenergy_eV[-1] * 1e3
        self.h_kenergy_rec_meV = self.h_kenergy_eV[-1] * 1e3
        self.total_kenergy_rec_meV = self.total_kenergy_eV[-1] * 1e3
        
        self.approx_sb = (self.total_kenergy_rec_meV + self.e_bind * 1e3 / q_e) / (self.thz_freq * 1e-12 / 0.2418)
        self.results_dict['individual total KE (meV)'] = '{:.2f}'.format(self.total_kenergy_rec_meV)
        self.results_dict['individual electron KE (meV)'] = '{:.2f}'.format(self.e_kenergy_rec_meV)
        self.results_dict['individual hole KE (meV)'] = '{:.2f}'.format(self.h_kenergy_rec_meV)
        self.results_dict['approximate sideband'] = '{:.1f}'.format(self.approx_sb)

        self.ion_time = (self.phi - np.pi / 2) / (2 * np.pi * self.thz_freq)
        self.recol_time = self.ion_time + self.recollision_time_ps * 1e-12


        self.results_dict['ionization time'] = '{:.3f}'.format(self.ion_time*1e15)
        self.results_dict['recollision time'] = '{:.3f}'.format(self.recol_time*1e15)

    def make_energy_slew_rate(self, time, splitting):
        """
        Ideally this will help calculate the Landau-Zener tunneling rate
        :param time:
        :return:
        """
        time = time * 1e-15 # now in seconds
        splitting = splitting * 1e-3 * q_e # Now in J
        h_kenergy_J_interp = spi.interp1d(self.time_vals_s, q_e * self.h_kenergy_eV)
        #print "Hole energy", h_kenergy_J_interp(time) / q_e
        #print "Hole k vector = {:.2e}".format(np.sqrt(2 * h_kenergy_J_interp(time) * self.m_h_eff) / hbar)
        dk_dt = -q_e * self.thz_field * np.cos(2 * np.pi * self.thz_freq * time + self.phi)
        #print "dk_dt = ", dk_dt
        dE_dk = np.sqrt(2 * h_kenergy_J_interp(time) / self.m_h_eff)  # This is now in units of velocity
        #print "dE_dk = ", dE_dk
        self.dE_dt_J_s = dE_dk * dk_dt # in SI
        self.energy_slew_rate_meV_fs = self.dE_dt_J_s * (1e3 / q_e) * 1e-15
        self.diabatic_probability = np.exp(-2 * np.pi * splitting**2 / (hbar * abs(self.dE_dt_J_s)))

    def make_parameters(self):
        """
        Creates the phase-independent attributes.
        """
        self.keldysh = (2 * np.pi * self.thz_freq) * np.sqrt(2 * self.reduced_mass * self.e_bind) / (q_e * self.thz_field)
        self.ponderomotive_J = (q_e**2 * self.thz_field**2) / (4 * self.reduced_mass * (2 * np.pi * self.thz_freq)**2)
        self.ponderomotive_meV = self.ponderomotive_J / q_e * 1e3 # Now its in meV
        self.max_sb_order = (3.17 * self.ponderomotive_J + self.e_bind) / (self.thz_freq * 6.626e-34)
        self.results_dict['keldysh parameter'] = '{:.3f}'.format(self.keldysh)
        self.results_dict['ponderomotive energy (meV)'] = '{:.2f}'.format(self.ponderomotive_meV)
        self.results_dict['maximum order sideband'] = '{:.1f}'.format(self.max_sb_order)

if __name__ == '__main__':
    from PyQt5 import QtGui
    app = QtGui.QApplication(sys.argv)
    ex = Window()
    #ex.show()
    sys.exit(app.exec_())
