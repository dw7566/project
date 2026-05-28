import xml.etree.ElementTree as ET
import pandas as pd
import os



def csv_mod(filename,custom_csv):

    tree = ET.parse(filename)

    # I at -1V[A]
    c = tree.find(".ElectroOpticalMeasurements/ModulatorSite/Modulator/PortCombo/IVMeasurement/Current")
    y_2 = c.text.split(",")
    y_list = list(map(float, y_2))
    y_list_1 = []
    for i in range(len(y_list)):
        g = abs(y_list[i])
        y_list_1.append(g)


    Lot = tree.find("./TestSiteInfo").get('Batch')
    Wafer = tree.find("./TestSiteInfo").get('Wafer')
    Mask = tree.find("./TestSiteInfo").get('Maskset')
    Column = tree.find("./TestSiteInfo").get('Diecolumn')
    Row = tree.find("./TestSiteInfo").get('DieRow')

    Name = tree.find('./ElectroOpticalMeasurements/ModulatorSite/Modulator/DeviceInfo').get("Name")

    df = pd.DataFrame(columns=['Lot', 'Wafer', 'Mask', 'Name','Row', 'Column', 'I at -1V [A]', 'I at 1V [A]'])

    df.loc[0] = [Lot, Wafer, Mask,  Name, Row, Column, y_list_1[4], y_list_1[12]]



    if custom_csv == 'T':
        location = 'process_Result.csv'
        # if datetime.datetime.now().minute in os.path.basename('.\\res\\csv\\'):
        if not os.path.exists('.\\res\\csv\\%s'%location):
            try:
                os.makedirs('.\\res\\csv\\')
            except:
                FileExistsError
                pass
            df.to_csv('.\\res\\csv\\%s'%location, mode='w', index=False)
        else:
            df.to_csv('.\\res\\csv\\%s'%location, mode='a', index=False, header=False)