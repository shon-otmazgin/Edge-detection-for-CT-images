classdef LineDisplay
    %LINEDISPLAY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        
        function lineStrctCell = calcLineStrctCellVariability(im,seg,prior,varMask)
            [~,overlay] = Utils.displayCertaintyUncertainty2_3D(im,prior,Utils.getBoundries(seg,false));
            lineStrctCell = LineStruct.getLineStructCell(overlay);
            lineStrctCell = LineStruct.getLineStructCell(varMask,[0,0,1],lineStrctCell);
        end
        
        function displayVariability(im,seg,prior,varMask)
            lineStrctCell = LineDisplay.calcLineStrctCellVariability(im,seg,prior,varMask);
            LineDisplay.ui(im,lineStrctCell);
        end
        
        
        function displaySegsOverlay(im,segs)
            if ~iscell(segs)
                temp = segs;
                segs = {temp};
            end
            
            colors = {[0,0,1],[0,0.5,1],[0,1,1],[0,1,0.75],[0.5,1,0]};
            if length(colors)<length(segs)
                error 'too many segs, not enough colors. function should be fixed';
            end
            
            lineStrctCell = LineStruct.getLineStructCell(segs{1},colors{1});
            for t=2:length(segs)
                lineStrctCell = LineStruct.getLineStructCell(segs{t},colors{2},lineStrctCell);
            end
            
            LineDisplay.ui(im,lineStrctCell);
        end
        
        function saveVariabilityAsTiff(im,seg,prior,varMask,tiffName)
            if exist(tiffName,'file')
                delete(tiffName);
            end
            
            firstSlice = find(sum(sum(seg,1),2) > 0,1,'first');
            nSlices = find(sum(sum(seg,1),2) > 0,1,'last') - firstSlice + 1;
            lineStrctCell = LineDisplay.calcLineStrctCellVariability(im,seg,prior,varMask);
            
            f = figure('units','normalized','outerposition',[0 0 1 1]);
            
            ax = axes('Units','pixels');
            set(ax, 'units', 'normalized', 'position', [0.05 0.15 0.9 0.8])
            for z = 1:nSlices
                [z, nSlices]
                
                LineDisplay.showZSliceData(im, lineStrctCell,z + firstSlice - 1,[min(im(:)),max(im(:))]);
                currentRes = getframe;
                
                if ~exist('outMat','var')
                    [m, n, ~] = size(currentRes.cdata);
                    outMat = uint8(zeros(m,n,nSlices,3));
                end
                outMat(:,:,z,:) = permute(currentRes.cdata,[1,2,4,3]);
            end
            close all;
            IO.saveAsTiff(tiffName,outMat);
        end
        
        function ui(im,lineStructCell)
            
            global g_im;
            global g_lineStrctCell;
            
            
            % Create a figure and axes
            f = figure('Visible','off','units','normalized','outerposition',[0 0 1 1]);
            ax = axes('Units','pixels');
            set(ax, 'units', 'normalized', 'position', [0.05 0.15 0.9 0.8])
            g_im = im;
            
            if exist('lineStructCell','var')
                g_lineStrctCell = lineStructCell;
            else
                g_lineStrctCell = [];
            end
            
            % Create pop-up menu
            %popup = uicontrol('Style', 'popup',...
            %       'String', {'parula','jet','hsv','hot','cool','gray'},...
            %       'Position', [20 340 100 50],...
            %       'Callback', @setmap);
            
            % Create push button
            %btn = uicontrol('Style', 'pushbutton', 'String', 'Clear',...
            %    'Position', [20 20 50 20],...
            %    'Callback', 'cla');
            
            % Create slider
            z = round(size(g_im,3)/2);
            sld = uicontrol('Style', 'slider',...
                'Min',1,'Max',size(g_im,3),'Value',z,...
                'Position', [500 50 120 20],...
                'Callback', @updateSliceCallback);
            
            % Add a text uicontrol to label the slider.
            txt = uicontrol('Style','text',...
                'Position',[500 75 120 20],...
                'String','Vertical Exaggeration');
            
            sldContrastMin = uicontrol('Style', 'slider',...
                'Min',min(g_im(:)),'Max',max(g_im(:)),'Value',min(g_im(:)),...
                'Position', [700 55 120 20],...
                'Callback', @updateContrastCallback);
            
            sldContrastMax = uicontrol('Style', 'slider',...
                'Min',min(g_im(:)),'Max',max(g_im(:)),'Value',max(g_im(:)),...
                'Position', [700 30 120 20],...
                'Callback', @updateContrastCallback);
            
            uicontrol('Style', 'pushbutton', 'String', 'contrast',...
                'Position', [900 55 120 20],...
                'Callback', @optContrastCallback);
            
            uicontrol('Style', 'pushbutton', 'String', 'moveChunk',...
                'Position', [900 30 120 20],...
                'Callback', @goToNextSegChunkCallback);
            
            set(gcf,'ResizeFcn', @figureResized)
            
            
            
            txtContrast = uicontrol('Style','text',...
                'Position',[700 80 120 20],...
                'String',getContrastTxt());
            
            % Make figure visble after adding all components
            f.Visible = 'on';
            % This code uses dot notation to set properties.
            % Dot notation runs in R2014b and later.
            % For R2014a and earlier: set(f,'Visible','on');
            
            
            minVal = round(sldContrastMin.get('Value'));
            maxVal = round(sldContrastMax.get('Value'));
            LineDisplay.showZSliceData(g_im, g_lineStrctCell,z,[minVal,maxVal]);
            txt.set('String',getTxt(size(g_im,3),z));
            %  function setmap(source,event)
            %      val = source.Value;
            %      maps = source.String;
            %      % For R2014a and earlier:
            %      % val = get(source,'Value');
            %      % maps = get(source,'String');
            
            %      newmap = maps{val};
            %      colormap(newmap);
            %  end
            
            function optContrastCallback(source,event)
                [~,range,~] = Utils.optimizeImContrast(g_im);
                sldContrastMin.set('Value',range(1))
                sldContrastMax.set('Value',range(2))
                updateContrastCallback(source,event);
            end
            
            function updateSliceCallback(source,event)
                z = round(sld.get('Value'));
                sld.set('Value',z);
                minVal = round(sldContrastMin.get('Value'));
                maxVal = round(sldContrastMax.get('Value'));
                txt.set('String',getTxt(size(g_im,3),z));
                LineDisplay.showZSliceData(g_im, g_lineStrctCell,z,[minVal,maxVal]);
                
            end
            
            
            
            function updateContrastCallback(source,event)
                
                minVal = round(sldContrastMin.get('Value'));
                maxVal = round(sldContrastMax.get('Value'));
                minVal = min(minVal,maxVal);
                maxVal = max(minVal,maxVal);
                sldContrastMin.set('Value',minVal);
                sldContrastMax.set('Value',maxVal);
                txtContrast.set('String',getContrastTxt());
                LineDisplay.showZSliceData(g_im, g_lineStrctCell,z, [minVal,maxVal]);
            end
            
            function str = getTxt(totalZ,z)
                str = [num2str(z) '/' num2str(totalZ)];
            end
            
            function [minVal,maxVal] = getMinMaxContastVals()
                minVal = round(sldContrastMin.get('Value'));
                maxVal = round(sldContrastMax.get('Value'));
            end
            function str = getContrastTxt()
                
                [minVal,maxVal] = getMinMaxContastVals();
                str = ['[' num2str(minVal) ',' num2str(maxVal) ']'];
            end
            
            function goToNextSegChunkCallback(object, eventdata)
                
                %extract line segments
                segMask = ~cellfun(@isempty,g_lineStrctCell);
                lineSegments = bwconncomp(~cellfun(@isempty,g_lineStrctCell));
                z  = sld.get('Value');
                
                %find current segment
                isZExists = @(arr)sum(arr==z)>0;
                if segMask(z)
                    currSeg = find(cellfun(isZExists,lineSegments.PixelIdxList),1,'first');
                    nextSeg = mod(currSeg,length(lineSegments.PixelIdxList))+1;
                    nextSlice = lineSegments.PixelIdxList{nextSeg}(1);
                else
                    nextSlice = z+find(segMask(z+1:end),1,'first');
                    if isempty(nextSlice)
                        nextSlice = find(segMask(1:z-1),1,'first');
                    end
                end
                sld.set('Value',nextSlice);
                updateSliceCallback();
            end
            % -=< Figure resize callback function >=-
            function figureResized(object, eventdata)
                FigPos = get(gcf,'Position');
                S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
                Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
                %   BtnStPnt = uint16(FigPos(3)-250)+1;
                %   if BtnStPnt < 300
                %       BtnStPnt = 300;
                %   end
                %   Btn_Pos = [BtnStPnt 20 100 20];
                %   ChBx_Pos = [BtnStPnt+110 20 100 20];
                %   if sno > 1
                %       set(shand,'Position', S_Pos);
                %   end
                %   set(stxthand,'Position', Stxt_Pos);
                %   set(ltxthand,'Position', Ltxt_Pos);
                %   set(wtxthand,'Position', Wtxt_Pos);
                %   set(lvalhand,'Position', Lval_Pos);
                %   set(wvalhand,'Position', Wval_Pos);
                %   set(Btnhand,'Position', Btn_Pos);
                %   set(ChBxhand,'Position', ChBx_Pos);
            end
        end
        
        
        
        function showZSliceData(im, lineStrctCell,z, rangeLim)
            imshow(im(:,:,z),rangeLim);
            if ~isempty(lineStrctCell) && ~isempty(lineStrctCell{z})
                lineStrctCell{z}.plotLines;
            end
        end
        
    end
    
end

