pro widgettest,dummy
;
; simple routine to test moving the base widget and extracting the new location.
;

common wids, base, mainstop

xlocation = 100
ylocation = 200

base = WIDGET_BASE(TITLE='Test Program for Widget Repositioning', $
                        /COLUMN, $
                        XSIZE=700, $
                        /ALIGN_LEFT, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='widgettest_event', $
                        XOFFSET=xlocation, $
                        YOFFSET=ylocation)

print,' base ID ',base


block1 = WIDGET_BASE(base, $
		/FRAME, $
		/COLUMN)

mainstop = WIDGET_BUTTON(block1, $
                         VALUE='Quit', $
                         UVALUE='QUIT', $
                         EVENT_PRO='widgettest_event', $
                         SENSITIVE=1, $
                         /FRAME)

print,' quit ID ',mainstop

WIDGET_CONTROL, base, /REALIZE

XMANAGER,"widgettest", base, /NO_BLOCK

end




pro widgettest_event, event

common wids, base, mainstop

if (event.id eq base) then begin
  print,' New coordinates = ',event.x,event.y
end else begin
  WIDGET_CONTROL, event.id, GET_UVALUE = eventval 

  print,' eventval = ',eventval

  CASE eventval OF
    'QUIT': 		WIDGET_CONTROL, base, /DESTROY
  ENDCASE
endelse

end

