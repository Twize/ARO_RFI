# ARO_RFI
Radio interference flagging/excision algorithms and associated functions for the Algonquin Radio Observatory (ARO) - Research project for UofT's AST425 (Supervisor: Prof. Keith Vanderlinde)

Files contained in this repository:
- dataIO.py (For loading in data from CITA Sunnyvale cluster)
- decimate.py (For decimating the raw data sets and outputting arrays of power and power^2) 
- radiometer.py (Quality metric for assessing effectiveness of RFI flagging/excision methods)
- upchannelize.py (Up-channelizing routine to remove narrowband persistent RFI prior to decimation) 
- SK.py (Spectral Kurtosis flagging algorithm)
- SumThreshold.py (SumThreshold iterative flagging algorithm, first developed in C++ by Andre Offringa) 

For any inquiries, contact me at: tyler.wizenberg@mail.utoronto.ca
