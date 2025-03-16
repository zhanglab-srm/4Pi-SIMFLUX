function saveFigure(figid, filepath)
% function saveFigure(figid, filepath)
% Save Figure as PNG file
    print(figid, '-dpng', '-r300', filepath);
end