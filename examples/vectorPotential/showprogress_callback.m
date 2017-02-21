function status=showprogress_callback(t,y,flag)

persistent h deltaT

status=0;

switch flag
    case 'init'
        h = waitbar(0,'1','Name','Solving system ...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0);
        deltaT = t(end)-t(1);
        

    case []
        waitbar(t/deltaT,h,sprintf('%f',t));
        if getappdata(h,'canceling')
            status=1;
        end

    case 'done'
        delete(h);
end
end
