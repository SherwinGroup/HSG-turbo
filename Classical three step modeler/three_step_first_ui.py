# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'three_step_first.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(916, 633)
        self.horizontalLayout = QtGui.QHBoxLayout(Form)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.frame = QtGui.QFrame(Form)
        self.frame.setMaximumSize(QtCore.QSize(200, 16777215))
        self.frame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtGui.QFrame.Raised)
        self.frame.setObjectName(_fromUtf8("frame"))
        self.gridLayout = QtGui.QGridLayout(self.frame)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_2 = QtGui.QGridLayout()
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.lineEdit_h_ke_rec = QtGui.QLineEdit(self.frame)
        self.lineEdit_h_ke_rec.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_h_ke_rec.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_h_ke_rec.setObjectName(_fromUtf8("lineEdit_h_ke_rec"))
        self.gridLayout_2.addWidget(self.lineEdit_h_ke_rec, 14, 1, 1, 1)
        self.label_e_ke_rec = QtGui.QLabel(self.frame)
        self.label_e_ke_rec.setObjectName(_fromUtf8("label_e_ke_rec"))
        self.gridLayout_2.addWidget(self.label_e_ke_rec, 13, 0, 1, 1)
        self.lineEdit_e_ke_rec = QtGui.QLineEdit(self.frame)
        self.lineEdit_e_ke_rec.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_e_ke_rec.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_e_ke_rec.setObjectName(_fromUtf8("lineEdit_e_ke_rec"))
        self.gridLayout_2.addWidget(self.lineEdit_e_ke_rec, 13, 1, 1, 1)
        self.label_h_ke_rec = QtGui.QLabel(self.frame)
        self.label_h_ke_rec.setObjectName(_fromUtf8("label_h_ke_rec"))
        self.gridLayout_2.addWidget(self.label_h_ke_rec, 14, 0, 1, 1)
        self.lineEdit_tot_ke_rec = QtGui.QLineEdit(self.frame)
        self.lineEdit_tot_ke_rec.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_tot_ke_rec.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_tot_ke_rec.setObjectName(_fromUtf8("lineEdit_tot_ke_rec"))
        self.gridLayout_2.addWidget(self.lineEdit_tot_ke_rec, 12, 1, 1, 1)
        self.label_total_ke_rec = QtGui.QLabel(self.frame)
        self.label_total_ke_rec.setObjectName(_fromUtf8("label_total_ke_rec"))
        self.gridLayout_2.addWidget(self.label_total_ke_rec, 12, 0, 1, 1)
        self.label_recollision_time = QtGui.QLabel(self.frame)
        self.label_recollision_time.setObjectName(_fromUtf8("label_recollision_time"))
        self.gridLayout_2.addWidget(self.label_recollision_time, 11, 0, 1, 1)
        self.lineEdit_recollision_time = QtGui.QLineEdit(self.frame)
        self.lineEdit_recollision_time.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_recollision_time.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_recollision_time.setObjectName(_fromUtf8("lineEdit_recollision_time"))
        self.gridLayout_2.addWidget(self.lineEdit_recollision_time, 11, 1, 1, 1)
        self.label_rec_group = QtGui.QLabel(self.frame)
        self.label_rec_group.setAlignment(QtCore.Qt.AlignCenter)
        self.label_rec_group.setObjectName(_fromUtf8("label_rec_group"))
        self.gridLayout_2.addWidget(self.label_rec_group, 9, 0, 1, 2)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem, 15, 0, 1, 1)
        self.label_m_e_eff = QtGui.QLabel(self.frame)
        self.label_m_e_eff.setMinimumSize(QtCore.QSize(100, 0))
        self.label_m_e_eff.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_m_e_eff.setObjectName(_fromUtf8("label_m_e_eff"))
        self.gridLayout_2.addWidget(self.label_m_e_eff, 0, 0, 1, 1)
        self.lineEdit_m_e_eff = QtGui.QLineEdit(self.frame)
        self.lineEdit_m_e_eff.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_m_e_eff.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_m_e_eff.setObjectName(_fromUtf8("lineEdit_m_e_eff"))
        self.gridLayout_2.addWidget(self.lineEdit_m_e_eff, 0, 1, 1, 1)
        self.label_keldysh = QtGui.QLabel(self.frame)
        self.label_keldysh.setMinimumSize(QtCore.QSize(100, 0))
        self.label_keldysh.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_keldysh.setObjectName(_fromUtf8("label_keldysh"))
        self.gridLayout_2.addWidget(self.label_keldysh, 16, 0, 1, 1)
        self.lineEdit_keldysh = QtGui.QLineEdit(self.frame)
        self.lineEdit_keldysh.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_keldysh.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_keldysh.setObjectName(_fromUtf8("lineEdit_keldysh"))
        self.gridLayout_2.addWidget(self.lineEdit_keldysh, 16, 1, 1, 1)
        self.label_ponderomotive = QtGui.QLabel(self.frame)
        self.label_ponderomotive.setMinimumSize(QtCore.QSize(100, 0))
        self.label_ponderomotive.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_ponderomotive.setObjectName(_fromUtf8("label_ponderomotive"))
        self.gridLayout_2.addWidget(self.label_ponderomotive, 17, 0, 1, 1)
        self.lineEdit_thz_field = QtGui.QLineEdit(self.frame)
        self.lineEdit_thz_field.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_thz_field.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_thz_field.setObjectName(_fromUtf8("lineEdit_thz_field"))
        self.gridLayout_2.addWidget(self.lineEdit_thz_field, 4, 1, 1, 1)
        self.lineEdit_thz_freq = QtGui.QLineEdit(self.frame)
        self.lineEdit_thz_freq.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_thz_freq.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_thz_freq.setObjectName(_fromUtf8("lineEdit_thz_freq"))
        self.gridLayout_2.addWidget(self.lineEdit_thz_freq, 3, 1, 1, 1)
        self.lineEdit_phase = QtGui.QLineEdit(self.frame)
        self.lineEdit_phase.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_phase.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_phase.setObjectName(_fromUtf8("lineEdit_phase"))
        self.gridLayout_2.addWidget(self.lineEdit_phase, 6, 1, 1, 1)
        self.label_thz_field = QtGui.QLabel(self.frame)
        self.label_thz_field.setMinimumSize(QtCore.QSize(100, 0))
        self.label_thz_field.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_thz_field.setObjectName(_fromUtf8("label_thz_field"))
        self.gridLayout_2.addWidget(self.label_thz_field, 4, 0, 1, 1)
        self.label_m_h_eff = QtGui.QLabel(self.frame)
        self.label_m_h_eff.setMinimumSize(QtCore.QSize(100, 0))
        self.label_m_h_eff.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_m_h_eff.setObjectName(_fromUtf8("label_m_h_eff"))
        self.gridLayout_2.addWidget(self.label_m_h_eff, 1, 0, 1, 1)
        self.lineEdit_max_sbs = QtGui.QLineEdit(self.frame)
        self.lineEdit_max_sbs.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_max_sbs.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_max_sbs.setObjectName(_fromUtf8("lineEdit_max_sbs"))
        self.gridLayout_2.addWidget(self.lineEdit_max_sbs, 18, 1, 1, 1)
        self.label_e_bind = QtGui.QLabel(self.frame)
        self.label_e_bind.setMinimumSize(QtCore.QSize(100, 0))
        self.label_e_bind.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_e_bind.setObjectName(_fromUtf8("label_e_bind"))
        self.gridLayout_2.addWidget(self.label_e_bind, 2, 0, 1, 1)
        self.label_thz_freq = QtGui.QLabel(self.frame)
        self.label_thz_freq.setMinimumSize(QtCore.QSize(100, 0))
        self.label_thz_freq.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_thz_freq.setObjectName(_fromUtf8("label_thz_freq"))
        self.gridLayout_2.addWidget(self.label_thz_freq, 3, 0, 1, 1)
        self.label_max_sbs = QtGui.QLabel(self.frame)
        self.label_max_sbs.setMinimumSize(QtCore.QSize(100, 0))
        self.label_max_sbs.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_max_sbs.setObjectName(_fromUtf8("label_max_sbs"))
        self.gridLayout_2.addWidget(self.label_max_sbs, 18, 0, 1, 1)
        self.lineEdit_e_bind = QtGui.QLineEdit(self.frame)
        self.lineEdit_e_bind.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_e_bind.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_e_bind.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lineEdit_e_bind.setObjectName(_fromUtf8("lineEdit_e_bind"))
        self.gridLayout_2.addWidget(self.lineEdit_e_bind, 2, 1, 1, 1)
        self.lineEdit_m_h_eff = QtGui.QLineEdit(self.frame)
        self.lineEdit_m_h_eff.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_m_h_eff.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_m_h_eff.setObjectName(_fromUtf8("lineEdit_m_h_eff"))
        self.gridLayout_2.addWidget(self.lineEdit_m_h_eff, 1, 1, 1, 1)
        self.label_arb_energy = QtGui.QLabel(self.frame)
        self.label_arb_energy.setMinimumSize(QtCore.QSize(100, 0))
        self.label_arb_energy.setObjectName(_fromUtf8("label_arb_energy"))
        self.gridLayout_2.addWidget(self.label_arb_energy, 7, 0, 1, 1)
        self.label_phi = QtGui.QLabel(self.frame)
        self.label_phi.setMinimumSize(QtCore.QSize(100, 0))
        self.label_phi.setMaximumSize(QtCore.QSize(100, 16777215))
        self.label_phi.setObjectName(_fromUtf8("label_phi"))
        self.gridLayout_2.addWidget(self.label_phi, 6, 0, 1, 1)
        self.lineEdit_energy_threshold = QtGui.QLineEdit(self.frame)
        self.lineEdit_energy_threshold.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_energy_threshold.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_energy_threshold.setObjectName(_fromUtf8("lineEdit_energy_threshold"))
        self.gridLayout_2.addWidget(self.lineEdit_energy_threshold, 7, 1, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(17, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem1, 5, 1, 1, 1)
        spacerItem2 = QtGui.QSpacerItem(17, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem2, 8, 0, 1, 1)
        self.lineEdit_ponderomotive = QtGui.QLineEdit(self.frame)
        self.lineEdit_ponderomotive.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_ponderomotive.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_ponderomotive.setObjectName(_fromUtf8("lineEdit_ponderomotive"))
        self.gridLayout_2.addWidget(self.lineEdit_ponderomotive, 17, 1, 1, 1)
        self.label_approx_sb = QtGui.QLabel(self.frame)
        self.label_approx_sb.setObjectName(_fromUtf8("label_approx_sb"))
        self.gridLayout_2.addWidget(self.label_approx_sb, 10, 0, 1, 1)
        self.lineEdit_approx_sb = QtGui.QLineEdit(self.frame)
        self.lineEdit_approx_sb.setMinimumSize(QtCore.QSize(50, 0))
        self.lineEdit_approx_sb.setMaximumSize(QtCore.QSize(55, 16777215))
        self.lineEdit_approx_sb.setObjectName(_fromUtf8("lineEdit_approx_sb"))
        self.gridLayout_2.addWidget(self.lineEdit_approx_sb, 10, 1, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_2, 0, 0, 1, 1)
        self.horizontalLayout.addWidget(self.frame)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.graphicsView_pos_graph = PlotWidget(Form)
        self.graphicsView_pos_graph.setObjectName(_fromUtf8("graphicsView_pos_graph"))
        self.verticalLayout_2.addWidget(self.graphicsView_pos_graph)
        self.graphicsView_energy_graph = PlotWidget(Form)
        self.graphicsView_energy_graph.setObjectName(_fromUtf8("graphicsView_energy_graph"))
        self.verticalLayout_2.addWidget(self.graphicsView_energy_graph)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem3)
        self.pushButton_save = QtGui.QPushButton(Form)
        self.pushButton_save.setObjectName(_fromUtf8("pushButton_save"))
        self.horizontalLayout_2.addWidget(self.pushButton_save)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.horizontalLayout.addLayout(self.verticalLayout_2)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)
        Form.setTabOrder(self.lineEdit_m_e_eff, self.lineEdit_m_h_eff)
        Form.setTabOrder(self.lineEdit_m_h_eff, self.lineEdit_e_bind)
        Form.setTabOrder(self.lineEdit_e_bind, self.lineEdit_thz_freq)
        Form.setTabOrder(self.lineEdit_thz_freq, self.lineEdit_thz_field)
        Form.setTabOrder(self.lineEdit_thz_field, self.lineEdit_phase)
        Form.setTabOrder(self.lineEdit_phase, self.lineEdit_energy_threshold)
        Form.setTabOrder(self.lineEdit_energy_threshold, self.lineEdit_approx_sb)
        Form.setTabOrder(self.lineEdit_approx_sb, self.lineEdit_recollision_time)
        Form.setTabOrder(self.lineEdit_recollision_time, self.lineEdit_tot_ke_rec)
        Form.setTabOrder(self.lineEdit_tot_ke_rec, self.lineEdit_e_ke_rec)
        Form.setTabOrder(self.lineEdit_e_ke_rec, self.lineEdit_h_ke_rec)
        Form.setTabOrder(self.lineEdit_h_ke_rec, self.lineEdit_keldysh)
        Form.setTabOrder(self.lineEdit_keldysh, self.lineEdit_ponderomotive)
        Form.setTabOrder(self.lineEdit_ponderomotive, self.lineEdit_max_sbs)
        Form.setTabOrder(self.lineEdit_max_sbs, self.pushButton_save)
        Form.setTabOrder(self.pushButton_save, self.graphicsView_pos_graph)
        Form.setTabOrder(self.graphicsView_pos_graph, self.graphicsView_energy_graph)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.label_e_ke_rec.setText(_translate("Form", "Electron KE (meV)", None))
        self.label_h_ke_rec.setText(_translate("Form", "Hole KE (meV)", None))
        self.label_total_ke_rec.setText(_translate("Form", "Total KE (meV)", None))
        self.label_recollision_time.setText(_translate("Form", "Time (ps)", None))
        self.label_rec_group.setText(_translate("Form", "At recollision:", None))
        self.label_m_e_eff.setText(_translate("Form", "m*_e", None))
        self.lineEdit_m_e_eff.setText(_translate("Form", "0.063", None))
        self.label_keldysh.setText(_translate("Form", "Keldysh", None))
        self.label_ponderomotive.setText(_translate("Form", "U_p (meV)", None))
        self.lineEdit_thz_field.setText(_translate("Form", "10", None))
        self.lineEdit_thz_freq.setText(_translate("Form", "0.54", None))
        self.lineEdit_phase.setText(_translate("Form", "60", None))
        self.label_thz_field.setText(_translate("Form", "F_THz (kV/cm)", None))
        self.label_m_h_eff.setText(_translate("Form", "m*_h", None))
        self.label_e_bind.setText(_translate("Form", "E_bind (meV)", None))
        self.label_thz_freq.setText(_translate("Form", "f_THz (THz)", None))
        self.label_max_sbs.setText(_translate("Form", "# SBs", None))
        self.lineEdit_e_bind.setText(_translate("Form", "10", None))
        self.lineEdit_m_h_eff.setText(_translate("Form", "0.51", None))
        self.label_arb_energy.setText(_translate("Form", "Energy line (meV)", None))
        self.label_phi.setText(_translate("Form", "Phase (deg)", None))
        self.lineEdit_energy_threshold.setText(_translate("Form", "10", None))
        self.label_approx_sb.setText(_translate("Form", "Approx. sideband", None))
        self.pushButton_save.setText(_translate("Form", "Save us!", None))

from pyqtgraph import PlotWidget
