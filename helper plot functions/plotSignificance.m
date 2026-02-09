function plotSignificance(stats,xPos,dataMax,scale)
% plot significance
sigValues = [0.05, 0.01, 0.001];

yLim = ylim;
yMax = yLim(end);
yRange = diff(yLim);
nPairs = size(stats.results,1);

hold on
for i = 1:nPairs
    pValue = stats.results{i,"P-value"};
    if pValue < sigValues(1)
        if nPairs == 1
            x1 = xPos(1);
            x2 = xPos(2);
        else
            switch cell2mat(stats.results{i,"Group A"})
                case cell2mat(stats.results{1,"Group A"})
                    x1 = xPos(1);
                case cell2mat(stats.results{1,"Group B"})
                    x1 = xPos(2);
                otherwise
                    x1 = xPos(3);
            end
    
            switch cell2mat(stats.results{i,"Group B"})
                case cell2mat(stats.results{1,"Group A"})
                    x2 = xPos(1);
                case cell2mat(stats.results{1,"Group B"})
                    x2 = xPos(2);
                otherwise
                    x2 = xPos(3);
            end   
        end
        text_x = x1 + 0.5*abs((x2-x1));

        switch scale
            case 'linear'
                horz_line_height_min = dataMax + 0.2*(yMax - dataMax);
                 if (abs(x2-x1) == 2)
                    horz_line_height = horz_line_height_min + 0.06*yRange;
                 else
                    horz_line_height = horz_line_height_min;
                 end
                 switch stats.sig
                     case 'exact'
                        text_y = horz_line_height + 0.04*yRange;
                     case 'symbol'
                        text_y = horz_line_height + 0.02*yRange;
                 end
            case 'log'
                horz_line_height_min = 10^( 0.2*( log10(yMax) - log10(dataMax) ) + log10(dataMax) );
                 if (abs(x2-x1) == 2)
                    horz_line_height = 10^( 0.06*( log10(yMax) - log10(yLim(1)) ) + log10(horz_line_height_min) );
                 else
                    horz_line_height = horz_line_height_min;
                 end
                switch stats.sig
                    case 'exact'
                       text_y = 10^( 0.04*( log10(yMax) - log10(yLim(1)) ) + log10(horz_line_height) );
                    case 'symbol'
                       text_y = 10^( 0.02*( log10(yMax) - log10(yLim(1)) ) + log10(horz_line_height) );
                end

        end

        plot( [x1 x2],horz_line_height*ones(1,2),'-k','Linewidth',1.2 )
        
        switch stats.sig 
            case 'exact'
                text_string = sprintf('p = %.3g',pValue);
                pvaluefontsize = 10;
            case 'symbol'
                if pValue < sigValues(1) && pValue > sigValues(2)
                    text_string = '*';
                elseif pValue < sigValues(2) && pValue > sigValues(3)
                    text_string = '**';
                else 
                    text_string = '***';
                end
                pvaluefontsize = 20;
        end

        text_box = text( text_x, text_y, text_string);
        text_box.FontSize = pvaluefontsize;
        text_box.FontWeight = 'bold';
        text_box.HorizontalAlignment = "center";
    end

end

hold off

end
