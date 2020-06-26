;This program automatically restarts the PIGI acquisition in auto_PIGI_underway_2020.vi

#include <AutoItConstants.au3>
#include <MsgBoxConstants.au3>

;Activate the LabVIEW (LV) windown that will run the restart program
   Local $hWnd = "Piggy_automation_restart"
   WinSetState($hWnd,"", @SW_SHOW)
   WinActivate("Piggy_automation_restart")

;Wait briefly
   Sleep (500)

;Message box - warn that the acquisition will be restarting soon
   MsgBox ($MB_ICONWARNING, " ", "Underway Acquisition restarting in 10 seconds", 10)

;Wait briefly
   Sleep (500)

;Activate LV window that runs the main acquisition
   Local $hWnd = "autoPIGI_underway_2020.vi"
   WinSetState($hWnd,"", @SW_SHOW)
   Sleep (1000)
   WinActivate("autoPIGI_underway_2020.vi")
   Sleep (1000)

;Adjust window size
   WinMove($hWnd,"",50,50,1100,700) ;position = left, bottom, width, height
   sleep(1000)

;Click LV stop button
   MouseMove(170, 120, 5) ;locate stop button
   sleep(1000)
   MouseClick($MOUSE_CLICK_LEFT)
   sleep(1000)

;Wait for LV acquisition to stop
   sleep(10000)

;Re-click the LV start button
   MouseMove(120, 120, 5) ;locate start button
   sleep(1000)
   MouseClick($MOUSE_CLICK_LEFT)

;Move the mouse away from the start arrow.
   MouseMove(150, 150, 5) ;

; Message box
   MsgBox ($MB_ICONWARNING, " ", "Please turn off the mouse when you are done with it", 5)
