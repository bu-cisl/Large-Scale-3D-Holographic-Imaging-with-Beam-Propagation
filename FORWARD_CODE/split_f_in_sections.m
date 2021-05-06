function f_split = split_f_in_sections(f,section_siz)
% used for split object into different sections for larege scale
% calculation

st = 1;
count = 1;
loop = 1;
[~,~,len] = size(f);
while loop
    en = st+section_siz-1;
    if en > len
        en = len;
        loop = 0;
    end
    f_split{count} = f(:,:,st:en);
    st = st + section_siz;
    count = count + 1;
end
end