%script to  read acquisition setup information
%from the xml file created using the acquisition software
%version >1.61
%
%uses the XMLparser toolbox
%http://www.artefact.tk/software/matlab/xml/
%
%input:
%       xml_filename:   xml filename to read
%
%output:
%
%   session_info:   structure storing the general information 
%                   about the experimental session.
%                   the structure contains the following fields:
%                       abstract: the abstract name
%                       protocol_code: the protocol code
%                       date: the date and hour of abstract creation
%                       place: the place where the session was conducted
%                       comments: general comments
%   sbj_info:       structure storing the patient information 
%                   the structure contains the following fields:
%                       id: patient unique identifier. ex: '24'
%                       family_name: ex: 'Uccelli'
%                       first_name: ex: 'Adriana'
%                       sex: 'F' or 'M'
%                       birth_date: ex: '1982-10-19 00:00:00'
%                       weight: ex: '50'
%                       height: ex: '168'
%                       registration_date: ex: '2006-11-10 00:00:00'
%                       pathology: 'none'
%                        comments: 'none'
%
%   board_info      structure storing the information 
%                   about the settings of the acquisition board
%                   the structure contains the following fields:
%                       fsamp:  sampling frequency (in Hz)
%                       ad_bits: number of bits of A/D converter
%
%   setup_info: structure storing the information 
%               about the setup used to acquire the signals
%               the structure contains the following fields:
%                  nchs: number of acquired channels. ex: 128
%                  ad_bits: number of bits of A/D converter
%                  narrays: number of arrays used in the setup
%                  array: vector with the characteristics of each array
%                         each element of this vector is a structure
%                         with the following fields:
%                            array_first_ch: first channel of the
%                               amplifier the array is connected to.
%                            muscle_code: code of the muscle where
%                               the array is placed. ex: 'UT'
%                            muscle_name: name of the muscle where
%                               the array is placed. ex: 'Upper trapezius'
%                            side: body side where the array is
%                               placed ('L', 'R' or '-')
%                            array: description of the array. ex:
%                               '12 el. Grid row - 8 mm IED'
%                            nchannels: number of channels used for
%                               the array
%
%       sig_info:       cell  array. One cell for each signal recorded 
%                       during the experimental session {sig}
%                       Each cell contains a structure
%                       The following fields are present:
%                       
%                           filename:   signal filename
%                           length:     length of the signal in seconds
%                           gain:       vector with one entry for each array
%                                       containing the gain used to acquire the
%                                       signal
%
%author: Marco Gazzoni
%date: 20070219
%
% Copyright (C) 27/10/2006 LISiN Politecnico di Torino (lisin@polito.it)
% All rights reserved

function [session_info, sbj_info, board_info, setup_info, sig_info] = load_xml(xml_filename)

cd xmltree
t= xmltree([xml_filename]);
s = convert(t);
cd ..

%obtain the information about the subject
session_info.abstract= s.abstract;
session_info.protocol_code= s.protocol_code;
session_info.date= s.time;
session_info.place= s.place;
session_info.comments= s.comments;

sbj_info= s.patient;
board_info.fsamp= str2num(s.board.fsample);
board_info.ad_bits= str2num(s.board.ad_bits);

%obtain the number of acquired signals

%obtain the number of acuired signals in the session
nsigs= length(s.signal);
if(nsigs==1)
    signals{1}= s.signal;
else
    signals= s.signal;
end

%read the setup
%number of sensors:
narrays= length(s.setup.ordering);
setup_info.nchs= str2num(s.board.channels);
setup_info.ad_bits= str2num(s.board.ad_bits);
setup_info.narrays= narrays;
index_ch= 1;
for array=1:narrays
    %obtain the first channel for the current array
    array_first_ch= str2num(s.setup.ordering{array}.number);
    setup_info.array(array).array_first_ch= array_first_ch;
    setup_info.array(array).muscle_code= s.setup.channel{index_ch}.muscle.code;
    setup_info.array(array).muscle_name= s.setup.channel{index_ch}.muscle.name;
    setup_info.array(array).side= s.setup.channel{index_ch}.side;
    setup_info.array(array).array= s.setup.channel{index_ch}.sensor.name;
    setup_info.array(array).nchannels= str2num(s.setup.channel{index_ch}.sensor.channels);
    index_ch= index_ch+setup_info.array(array).nchannels;
end

%cycle on the signals
for sig=1:nsigs
    %load the signal
    sig_info{sig}.filename= signals{sig}.file_name;
    sig_info{sig}.length= str2num(signals{sig}.duration);
    
    %obtain the gains used for the different detection systems
    for array=1:narrays
        g= signals{sig}.gain(array);
        sig_info{sig}.gain(array)= str2num(g{1});
    end
end




