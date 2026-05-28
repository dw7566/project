![header](https://capsule-render.vercel.app/api?type=wave&color=auto&height=200&section=header&text=SBBFU%20PROJECT&fontSize=50)

###### PE02_A02
### :pencil:Contents

[1. Introduction](#1-introduction)\
[2. Project information](#2-project-information)\
[3. Install and Run](#3-install-and-run)\
[4. Description of the module file feature](#4-description-of-the-module-file-feature)

---

# SBBFU Project :
##### Hi :wave:
##### Thank you for looking for us. The project is called SBBFU, which means "Special Black Box For U."
##### Analyze customer-provided data and provide results.

---

## 1. Introduction
#### - PE2 : SBBFU Project
We aim to develop Python-based data analysis software. \
The goal is to develop software that can analyze wafer scale data in detail.

#### - contributors : If you have any questions, please contact us at the following email.

|name|E-mail|:smile:|
|:---:|:---:|:---:|
|Park seoung min:boy:|psm401@hanyang.ac.kr|<img src = "https://user-images.githubusercontent.com/84078034/121825519-e074d280-cced-11eb-8507-995dd571eaf4.png" width = "100" height = "100">|
|Kim kang seok:boy:|stone410@naver.com|<img src = "https://user-images.githubusercontent.com/89181329/132611009-a175bf6e-38ad-4a77-844d-20c256de744b.jpg" width = "100" height = "100">|

---

## 2. Project information
 #### 
 + **Detailed project**
 
     The main task is to receive data from the customer and then receive the request.
     Main request 
     + Wafer
     + Die row & colum
     + option
 
   We can analyze the selected data by specifying specific elements
   In order to develop such analysis software, a black box(software) that implements this function is created to solve the customer's request.

 #### 
+ **Run file description**

   First, extract only the file named 'LMZ' from the file that the customer gave us.\
   Then, after analyzing the raw data given by the customer using various modules, make a csv file.

---

## 3. Install and Run

####
* Getting Stared
   + Entered the Terminal, write down 'pip install -r requirements.txt' and download it. \
``(base) C:\Download\PE02_A02_SBBFU-1>pip install -r requirements.txt``

* How to Run
  + Choose the raw data folder customer want to analyze.
   
  + After select the data you want to analyze, select various options such as wafer, die option, figure shape (show figure, save figure, save csv) and press run button.

---

## 4. Description of the module file feature

* Fitting module
   + The graph is drawn by parsing the raw data of Wavelengthsweep, IL, Current, and Voltage in the xml file.
   + The fitting of parsing a raw data and displays the data value y-axis corresponding to x-axis by the customer and the desired R-squared, etc. and stored in the graph to visualize the image.

* CSV module
  + It contains a variety of data information, including Lot, Wafer, and Operator etc.
  + Create a dataframe so that the meausured information in the xml file can be viewed at a glance.
  + Save this data frome in csv format in the 'Result' folder.
 
 ---
## 5. Run file algorithm
* First of all 
   + Please select options such as Data location, Analyze folder, Figure shape, Colum & Row, Save figure mode, Show figure mode, Save CSV in the run.py. 
   

* And then
  + The 'run' function in the runfile.py is executed.
* Next
   + Depending on the option initially selected, ivfitting.py, measured.py, spectra.py, process.py, fitting.py are executed. 
 ---

![image](https://user-images.githubusercontent.com/77437180/132507266-4615ae28-653f-494b-85aa-2b714852021e.png)
 
### :warning:precautions

 1) You must select all the options in the 'Checkbox'. If you don't, you'll make a error.
 2) When fitting sometimes doesn't work, it keeps running until fitting is done well.
 3) If you choose the All option when Wafer is selected, you must not choose any other option.
 
