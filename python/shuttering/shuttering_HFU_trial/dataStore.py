# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 11:13:54 2019

@author: Sarah
"""
class Datasets():

    def __init__(self, DATASET):

        if(DATASET == 0):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V00/'
            self.files = [
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_11-05-14_Trial-1.csv',
                   # 'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_11-07-40_Trial-1.csv',
                   # 'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_11-09-23_Trial-1.csv',
                   # 'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_11-11-30_Trial-1.csv',
                     ]
                # For data
            self.xLabels = ['Normal','12.5mm','10.5mm', '9.5mm']




        if(DATASET == 1):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V01/'
            self.files = [
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-24-11_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-28-35_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-31-07_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-34-46_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-37-21_Trial-1.csv',
                    #'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-40-34_Trial-1.csv', # electrode fell off and laugh/cough
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-43-07_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-45-28_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-47-18_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-54-06_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-55-50_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_01-58-30_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_02-00-32_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_02-02-40_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_02-04-34_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_02-07-33_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-01_02-09-41_Trial-1.csv',
                     ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']




        if(DATASET == 2):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V02/'
            self.files = [
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-10-10_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-16-12_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-19-34_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-22-01_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-24-55_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-27-18_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-29-42_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-31-54_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-38-06_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-40-29_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-43-32_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-46-29_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-49-45_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-52-07_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-55-23_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-02_10-57-25_Trial-1.csv',
                     ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']




        if(DATASET == 3):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V03/'
            self.files = [
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-33-17_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-36-29_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-39-26_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-41-47_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-44-17_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-46-51_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-50-39_Trial-1.csv',
#                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-53-19_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-57-50_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_08-59-34_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-01-36_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-03-33_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-05-28_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-05-56_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-09-54_Trial-1.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-11-51_Trial-1.csv',
                     ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']




        if(DATASET == 4):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V04/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-20-47_Trial-1.csv', # Laugh
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-22-28_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-24-10_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-26-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-27-40_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-29-46_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-31-45_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-33-43_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-35-39_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-40-51_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-42-50_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-44-53_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-46-30_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-48-46_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-50-47_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-53-36_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-04_09-55-53_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']



        if(DATASET == 5):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V05/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-32-17_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-35-20_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-38-17_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-40-47_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-44-09_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-46-00_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-48-07_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-49-49_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-54-34_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-56-27_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_11-58-19_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_12-00-32_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_12-02-51_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_12-05-29_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_12-07-25_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-08_12-09-05_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']


        if(DATASET == 6):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V06/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-02-42_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-06-49_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-17-59_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-19-43_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-21-11_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-23-09_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-24-50_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-26-37_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-28-07_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-33-16_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-35-27_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-37-29_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-39-16_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-41-17_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-43-07_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-45-14_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_01-47-08_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']


        if(DATASET == 7):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V07/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-02-16_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-04-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-06-27_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-08-17_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-10-23_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-12-31_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-14-53_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-17-00_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-21-56_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-23-42_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-25-55_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-28-18_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-30-46_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-33-02_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-35-49_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-09_02-37-55_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal', 'Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']


        if(DATASET == 8):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V08/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-08-36_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-10-32_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-12-40_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-14-23_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-16-42_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-18-40_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-20-58_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-22-51_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-28-31_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-30-54_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-31-39_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-36-31_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-38-56_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-41-07_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-43-47_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-46-04_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']


        if(DATASET == 9):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V09/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-52-26_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-54-50_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-56-55_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_11-58-32_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-00-13_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-01-44_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-03-15_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-04-47_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-06-24_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-07-40_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-09-13_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-11-05_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-13-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-14-39_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-16-29_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-18-04_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','Normal','Normal','12.5mm','12.5mm','12.5mm','12.5mm','10.5mm', '10.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','9.5mm','9.5mm']
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']


        if(DATASET == 10):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V10/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-23-32_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-25-09_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-27-02_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-28-31_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-30-10_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-31-41_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-33-18_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-34-45_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-37-34_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-39-19_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-41-12_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-43-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-44-58_Trial-1.csv',# cheek hold
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-47-43_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-49-56_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-51-49_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_12-53-14_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','Normal','Normal','12.5mm','12.5mm','12.5mm','12.5mm','10.5mm','10.5mm','10.5mm','10.5mm','10.5mm','9.5mm','9.5mm','9.5mm','9.5mm']
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']


        if(DATASET == 11):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V11/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-07-14_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-09-20_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-10-52_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-12-32_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-14-31_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-16-01_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-17-30_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-18-59_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-20-52_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-22-37_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-24-07_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-25-51_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-27-38_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-29-08_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-30-52_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-11_01-32-30_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','Normal','Normal','12.5mm','12.5mm','12.5mm','12.5mm','10.5mm', '10.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','9.5mm','9.5mm']
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']

        if(DATASET == 12):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V12/'
            self.files = [
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_10-50-58_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_10-52-36_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_10-54-39_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_10-55-46_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_10-56-59_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_10-58-15_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-00-15_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-01-47_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-06-21_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-07-44_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-09-04_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-10-33_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-11-58_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-13-24_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-15-14_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-14_11-16-42_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']


        if(DATASET == 13):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V13/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-19-52_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-22-15_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-24-09_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-25-40_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-27-36_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-29-17_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-30-55_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-32-21_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-37-55_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-39-21_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-40-59_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-42-36_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-44-23_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-45-57_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-47-59_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_11-49-39_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']

        if(DATASET == 14):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V14/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-01-16_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-03-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-05-06_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-06-41_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-08-50_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-10-27_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-12-40_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-14-26_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-19-51_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-21-23_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-23-24_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-25-09_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-27-23_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-29-03_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-31-01_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-15_12-33-03_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']


#        if(DATASET == 15):
#            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
#            self.path = './data/HFU-VTrials/V15/'
#            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-48-20_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-51-10_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-53-18_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-55-47_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-58-39_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-00-07_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-01-32_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-02-53_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-05-36_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-07-33_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-09-11_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-10-48_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-13-19_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-14-53_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-17-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-19-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-24-12_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-25-36_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-26-55_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-28-32_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-30-36_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-31-44_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-33-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-34-23_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-36-21_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-37-58_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-39-56_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-41-32_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-43-12_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-44-27_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-45-42_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-46-54_Trial-1.csv',
#                 ]
#            # For data
#            self.xLabels = ['Normal','Normal','Normal','Normal','12.5mm','12.5mm','12.5mm','12.5mm','10.5mm', '10.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','9.5mm','9.5mm','Normal','Normal','Normal','Normal','12.5mm','12.5mm','12.5mm','12.5mm','10.5mm', '10.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','9.5mm','9.5mm']



        if(DATASET == 15): # NO cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V15/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-48-20_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-51-10_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-58-39_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-00-07_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-05-36_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-07-33_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-13-19_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-14-53_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-24-12_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-25-36_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-30-36_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-31-44_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-36-21_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-37-58_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-43-12_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-44-27_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']




        if(DATASET == 150): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V15/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-53-18_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_02-55-47_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-01-32_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-02-53_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-09-11_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-10-48_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-17-08_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-19-08_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-26-55_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-28-32_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-33-08_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-34-23_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-39-56_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-41-32_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-45-42_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-10-30_03-46-54_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']



        if(DATASET == 16): # NO cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V16/'
            self.files = [
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-12-16_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-13-53_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-17-23_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-18-31_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-22-19_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-23-28_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-27-15_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-28-22_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-34-56_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-36-09_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-39-44_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-41-00_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-44-25_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-45-23_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-49-16_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-50-21_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']




        if(DATASET == 160): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V16/'
            self.files = [
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-15-02_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-16-07_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-19-48_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-20-55_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-24-40_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-25-43_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-29-32_Trial-1.csv',
#'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-30-38_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-37-14_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-38-15_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-41-59_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-42-58_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-46-34_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-47-38_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-51-25_Trial-1.csv',
'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_02-52-29_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']



        if(DATASET == 17): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V17/'
            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-03-11_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-04-35_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-09-25_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-10-39_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-15-09_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-16-28_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-22-26_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-23-39_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-30-58_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-32-25_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-36-44_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-37-10_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-42-32_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-43-01_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-48-21_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-49-42_Trial-1.csv',








                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']



        if(DATASET == 170): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V17/'

            self.files = [
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-06-26_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-07-45_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-11-54_Trial-1.csv',
#                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-13-10_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-17-55_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-19-23_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-25-01_Trial-1.csv',
 #               'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-26-41_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-33-56_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-35-12_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-39-36_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-40-52_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-45-19_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-46-39_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-52-32_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-06_03-53-45_Trial-1.csv',

                 ]
            # For data
            self.xLabels = ['Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm','Normal','Normal','12.5mm','12.5mm','10.5mm', '10.5mm','9.5mm','9.5mm']



        if(DATASET == 18): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V18/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-29_01-39-03_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']

        if(DATASET == 180): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V18/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-11-29_01-42-21_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']

        if(DATASET == 112): # WITH cheeckholding
            self.path = './data/HFU-VTrials/V11_cheekAndNot/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-02_10-41-56_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-02_10-43-38_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']

        if(DATASET == 1120): # WITH cheeckholding
            self.path = './data/HFU-VTrials/V11_cheekAndNot/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-02_10-45-08_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-02_10-46-36_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']



        if(DATASET == 19): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V19/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-09_11-29-41_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-09_11-32-05_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']

        if(DATASET == 190): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V19/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-09_11-34-38_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-09_11-37-05_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']


        if(DATASET == 20): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V20/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-10_10-15-17_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-10_10-16-59_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']

        if(DATASET == 200): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V20/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-10_10-19-06_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-10_10-20-54_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']


        if(DATASET == 21): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V21/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-10_01-49-06_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-10_01-51-31_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']

        if(DATASET == 210): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V21/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-10_01-54-03_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-10_01-56-25_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']


        if(DATASET == 22): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V22/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-11_02-22-49_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-11_02-24-18_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']

        if(DATASET == 220): # WITH cheeckholding
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU-VTrials/V22/'

            self.files = [
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-11_02-26-00_Trial-1.csv',
                'Corrected_data_RespiratoryOcclusionMeasurement_19-12-11_02-27-17_Trial-1.csv',
                 ]
            # For data
            self.xLabels = ['Normal','Normal']



        if(DATASET == 80):
            #path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.path = './data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
            self.files = [
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-06-19_11-15-23_Trial-1_NO.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-06-19_11-19-05_Trial-1_12.5mm.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-06-19_11-22-47_Trial-1_10.5mm.csv',
                    'Corrected_data_RespiratoryOcclusionMeasurement_19-06-19_11-26-04_Trial-1_9.5mm.csv',
                    #'Corrected_data_SlowSpirometryMeasurement_19-06-19_11-11-24_Trial-1.csv',
                     ]
                # For data
            self.xLabels = ['Normal','12.5mm','10.5mm', '9.5mm']


