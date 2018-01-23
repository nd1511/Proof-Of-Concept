
% this script tests plot_Fig4.p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by:           Nikolaos Dionelis
% Date created:         23 January 2018
% Project:              Phase-Aware Single-Channel Speech Enhancement using Modulation-Domain Kalman Filtering
% Supervisor:           M. Brookes
% Short Description:    Script to test plot_Fig4.p
%
% References:
%
%   [1] N. Dionelis and M. Brookes, "Phase-Aware
%       Single-Channel Speech Enhancement using Modulation-Domain
%       Kalman Filtering," Submitted to IEEE Trans. on Audio, Speech and
%       Language Process., 2017.
%
%   [2] M. Brookes, "VOICEBOX: A speech processing toolbox
%       for MATLAB," http://www.ee.ic.ac.uk/hp/staff/dmb/
%       voicebox/voicebox.html, 1997-2018.
%
% Calls:
% ======
%   plot_Fig4.p
%   (VOICEBOX)
%
% Inputs:
% =======
%   None
%
% Outputs:
% ========
%   None
%
% Instructions for use:
% =====================
%   Download: VOICEBOX from http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%   Download: babble.wav from http://www.speech.cs.cmu.edu/comp.speech/Section1/Data/noisex.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE WORK (AS DEFINED BELOW) IS PROVIDED UNDER THE TERMS OF CREATIVE
% COMMONS BY-ATTRIBUTION NON-COMMERCIAL LICENCE 3.0 WHICH CAN BE FOUND AT
% http://creativecommons.org/licenses/by-nc/3.0/legalcode.
% THE WORK IS PROTECTED BY COPYRIGHT AND/OR OTHER APPLICABLE LAW. ANY USE
% OF THE WORK OTHER THAN AS AUTHORIZED UNDER THIS LICENSE OR COPYRIGHT LAW
% IS PROHIBITED.
%
% BY EXERCISING ANY RIGHTS TO THE WORK PROVIDED HERE, YOU ACCEPT AND AGREE
% TO BE BOUND BY THE TERMS OF THIS LICENSE. TO THE EXTENT THIS LICENSE MAY
% BE CONSIDERED TO BE A CONTRACT, THE LICENSOR GRANTS YOU THE RIGHTS
% CONTAINED HERE IN CONSIDERATION OF YOUR ACCEPTANCE OF SUCH TERMS AND
% CONDITIONS.
%
% --------------------------------------------------------------------------------
% In plain English, under this licence you are free:
%
% To Share: to copy, distribute and transmit the work.
% to Remix: to adapt the work under the following conditions:
%
% - Attribution: You must attribute the work in the manner specified by the
% author or licensor (but not in any way that suggests that they endorse
% you or your use of the work).
% - Noncommercial: You may not use this work for commercial purposes.
%
% With the understanding that:
%
% Waiver: Any of the above conditions can be waived if you get permission
% from the copyright holder.
% Public Domain: Where the work or any of its elements is in the public
% domain under applicable law, that status is in no way affected by the
% license.
% Other Rights: In no way are any of the following rights affected by the
% license: your fair dealing or fair use rights, or other applicable
% copyright exceptions and limitations; the author's moral rights; rights
% other persons may have either in the work itself or in how the work is
% used, such as publicity or privacy rights.
% Notice: For any reuse or distribution, you must make clear to others the
% license terms of this work. The best way to do this is with a link to
% the URI above.
%
% -------------------------------------------------------------------------------
%
% Representations, Warranties and Disclaimer
%
% UNLESS OTHERWISE MUTUALLY AGREED TO BY THE PARTIES IN WRITING, LICENSOR
% OFFERS THE WORK AS-IS AND MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY
% KIND CONCERNING THE WORK, EXPRESS, IMPLIED, STATUTORY OR OTHERWISE,
% INCLUDING, WITHOUT LIMITATION, WARRANTIES OF TITLE, MERCHANTIBILITY,
% FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF
% LATENT OR OTHER DEFECTS, ACCURACY, OR THE PRESENCE OF ABSENCE OF ERRORS,
% WHETHER OR NOT DISCOVERABLE. SOME JURISDICTIONS DO NOT ALLOW THE
% EXCLUSION OF IMPLIED WARRANTIES, SO SUCH EXCLUSION MAY NOT APPLY TO YOU.
%
% Limitation on Liability.
%
% EXCEPT TO THE EXTENT REQUIRED BY APPLICABLE LAW, IN NO EVENT WILL
% LICENSOR BE LIABLE TO YOU ON ANY LEGAL THEORY FOR ANY SPECIAL,
% INCIDENTAL, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES ARISING OUT OF
% THIS LICENSE OR THE USE OF THE WORK, EVEN IF LICENSOR HAS BEEN ADVISED
% OF THE POSSIBILITY OF SUCH DAMAGES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add the datapath to voicebox
addpath ./voicebox

% Add the datapath to the external functions
addpath ./external_functions

%for snr_value = 15 : 5 : 15
for snr_value = 10 : 5 : 10
    
    snr_value
    
    for filename2 = {'SA1.WAV'}
    %for filename2 = {'SA2.WAV'}
        
        filename2
        
        % Read the PHN TIMIT file
        [main_true_reference_signal, fs] = readsph(filename2{1},'wt');
        
        timet=((1:length(main_true_reference_signal))-1)*(1/fs);
        
        %timet = timet(timet<=2);
        timet = timet(timet<=1);
        
        main_true_reference_signal = [zeros(length(timet),1); main_true_reference_signal];
        
        % Inner for loop
        for filename = {'white.wav'}
            
            filename
            
            % Read the noise file
            [s2, fs2] = readwav(filename{1});
            
            % Define the SNR
            snr = snr_value;
            
            s = v_addnoise(main_true_reference_signal,fs,snr,'',s2,fs2);
            
            % s is the noisy speech signal
            % main_true_reference_signal is the clean signal
            
            plot_Fig4(s);
            
        end
    end
end
