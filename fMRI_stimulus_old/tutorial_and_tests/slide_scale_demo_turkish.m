ListenChar(2);
% krem rengi bir pencere olu?turdum
[myWindow, rect] = Screen('OpenWindow',0,[230 230 230]);

%imleci saklad?m
Screen('HideCursorHelper', 0, 0)

%ekran?n merkezini belirledim
centerX=rect(3)/2;
centerY=rect(4)/2;

% imleci ekran?n mekerkezine koydum
SetMouse(centerX, centerY);

% while loop kullanarak ekrandaki di?er ögeleri olu?trum
buttons=0;
while ~buttons
    % mouse hareketlerini al?yorum
    [x, y, buttons] = GetMouse(myWindow);


    % x imlecin yatay eksendeki yerini belitiyor
    % 383 ve 983 ye?il dörtgenimin s?n?rlar? imlecimin yani k?rm?z?
    % dörtgenimin bu s?n?rlar içinde kalmas? için if kullanarak bir
    % statement olu?turuyorum
    if ( x > 383) & (x < 983)

    % ekrana ilgili cümleyi yazd?r?yorum    
    Screen('TextSize', myWindow,20);
    Screen('TextFont', myWindow, 'Times');
    [normBoundsRect, offsetBoundsRect]=Screen('TextBounds',myWindow, 'Adjust the slide below according to your confidence level');
    Screen('DrawText', myWindow, 'Adjust the slide below according to your confidence level', (centerX-(normBoundsRect(3)/2)),(centerY-(normBoundsRect(4)/2+100)), [0,0,0]);


    % ekran?n bir taraf?na not confident yazd?r?yorum
    Screen('TextSize', myWindow,20);
    Screen('TextFont', myWindow, 'Times');
    [normBoundsRect, offsetBoundsRect]=Screen('TextBounds',myWindow, 'Not Confident');
    Screen('DrawText', myWindow, 'Not Confident', (centerX-(normBoundsRect(3)/2-250)),(centerY-(normBoundsRect(4)/2+30)), [0,0,0]);


    % ekran?n di?er taraf?na very confident yazd?r?yorum
    Screen('TextSize', myWindow,20);
    Screen('TextFont', myWindow, 'Times');
    [normBoundsRect, offsetBoundsRect]=Screen('TextBounds',myWindow, 'Very Confident');
    Screen('DrawText', myWindow, 'Very Confident', (centerX-(normBoundsRect(3)/2+250)),(centerY-(normBoundsRect(4)/2+30)), [0,0,0]);


    % ye?il dörtgenimi çiziyorum, asl?nda belirledi?im centerlar?
    % kullanarak yazmay? denedim ancak çal??mad? o yüzden kendi ekran
    % boyutuma göre orta noktaya yerle?tirdim(1366*768)
    % Screen('FrameRect',myWindow,[50,205,50], [centerX-300, centerY +10,centerX+300, centerY-10],... 5);
     Screen('FrameRect',myWindow,[50,205,50], [383, 374,983, 394], [5])

     % ye?il dörtgenin içini beyaza boayad?m
     Screen('FillRect', myWindow,[255,255,255], [388,379,979,389])


    % ye?il kutunun alt?nda yazmas? gerekeni yazd?rd?m
    Screen('TextSize', myWindow,15);
    Screen('TextFont', myWindow, 'Times');
    [normBoundsRect, offsetBoundsRect]=Screen('TextBounds',myWindow, 'Note: Slide red rectangle to indicate your confidence level');
    Screen('DrawText', myWindow, 'Note: Slide red rectangle to indicate your confidence level', (centerX-(normBoundsRect(3)/2)),(centerY-(normBoundsRect(4)/2-100)), [0,0,0]);


     % x'i yani imlecin koordinatlar?n? kullanarak, 10 piksel kal?nl???nda,
     % x'in sa??na ve soluna 5 piksel gidecek ?ekilde k?rm?z? kutuyu
     % çizdirdim
    Screen('FillRect',myWindow,[255 50 50],[x-5,374, x+5, 394])

            % kutuyu ekranda gösterdim
            Screen('Flip',myWindow)
    end

    % yeni bir if kullanarak e?er farenin tu?lar?n?n bas?l?rsa bu loop'u
    % bitirmesini söyledim
    if buttons==1
        break
    end

end 


% clickleri toplamas?n? söyledim
[clicks, x, y, buttons] = GetClicks(myWindow);


% ListenChar kodunu kulland?m ama san?r?m i?e yaram?yor, ctrl+c yapmadan klavye
% kullan?lam?yor
ListenChar(1);


Screen('CloseAll');

% command window'a nereye t?kland???n? yazd?rd?m
fprintf('You clicked at %d and %d x-y location.\n', x, y);