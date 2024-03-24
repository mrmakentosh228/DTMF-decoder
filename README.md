# DTMF-decoder
### Краткое руководство по компиляции программы (язык C).
Программа была написана на C, используются такие библиотеки как: stdio, math
и stdlib. Для сборки проекта был использован CMake, для компиляции
программа CLion от JetBrains, стандарт языка C99. Для расшифровки
текстового файла с сигналом необходимо поместить интересующий файл в
одну директорию с файлом программы. Далее, запустить программу, ввести
ПУТЬ до интересующего файла с указанием расширения (.txt), будет выведена
декодированная последовательность символов.
### Краткое руководство по компиляции программы (python).
Программа была написана на python 3.8, необходимо наличие таких библиотек
как numpy и matplotlib (модуля pyplot), а также math. Для расшифровки
текстового файла с сигналом необходимо поместить интересующий файл в
одну директорию с файлом программы. Далее, запустить программу, ввести
название интересующего файла с указанием расширения (.txt), будет выведена
декодированная последовательность символов и график спектральной
амплитуды для восьми частот, используемых при кодировке сигнала.
### Что вообще программа делает, как работает.
Программа расшифровывает DTMF-сигнал используя следующую таблицу
(взято из википедии):

![изображение](https://github.com/mrmakentosh228/DTMF-decoder/assets/44507420/6b855ea2-0d97-47b4-a57b-cde3901bc253)

Считается дискретное преобразование фурье для каждого участка сигнала
(символа) с помощью [алгоритма Герцеля](https://ru.wikipedia.org/wiki/%D0%90%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC_%D0%93%D1%91%D1%80%D1%86%D0%B5%D0%BB%D1%8F) для каждой из указанных в таблице
частот, на основе полученных спектральных амплитуд делается вывод о том,
какой именно символ зашифрован на данном участке сигнала.
#### Краткое описание алгоритма.
Суть такова: есть два “ползущих” по сигналу буфера (в реализации на C это
очереди, на python - списки), один буфер для обнаружения тишины и пауз в 
сигнале, один для собственно вычисления ДПФ, он заполняется теми отсчетами 
(значения уровня звукового сигнала в декодируемом файле), которые кодируют символ.
Перебирая каждый элемент сигнала, сначала заполняется буфер паузы, если
обнаруживается пауза, то суммируется информация о всех найденных
символах в сигнале до этого момента, и на основе этого в конечную выходную
(декодированную последовательность) добавляется найденный символ. Пауза
обнаруживается на основе величины стандартного отклонения набора отсчетов
в буфере паузы (для тишины и переходного участка она меньше, чем для
участка с “полезным” сигналом). Если пауза не обнаружена, то отсчеты из
буфера паузы переходят в буфер для обнаружения символа. Как только
заполняется буфер символа, происходит расчет амплитуд для восьми частот
для участка сигнала, находящегося в буфере, на основе данной информации
вычисляется закодированный в данном участке символ, для данного символа
счетчик увеличивается на единицу, самые старые отсчеты в буфере
освобождаются для добавления новых. На участке сигнала, который кодирует
определенный символ, “символьный” буфер может заполняться несколько раз,
именно для этого считается кол-во появлений (счетчик) определенного символа
в результате вычислений на данном участке (кодирующем символ) сигнала.
При обнаружении паузы (или переходного участка) на основе данных счетчиков
и делается вывод о том, какой символ на данном участке закодирован.
#### Отличие реализации на C от реализации на python.
В силу отсутствия в языке си готовой реализации списков, было принято
решение использовать для буферов структуру очереди (первый вошел первый
вышел). Я думал об использовании динамических массивов для буферов, но
все же решил остановиться на очередях. В силу отсутствия словарей, для
определения символов использовались их индексы в соответствующих
массивах для частот. Также была дополнительно написана функция для поиска
индекса максимального элемента массива. Вследствие этих фактов в
реализации на си больше переменных по сравнению с реализацией на питоне.
#### Требования к сигналу.
В программе присутствуют числовые константы, которые получены методом
научного тыка, а именно длина буфера паузы и символа в секундах (0.007 сек и
0.08 сек соответственно), максимальная величина стандартного отклонения
для буфера паузы (0.32), минимальный порог для обнаружения спектральной
амплитуды определенной частоты (3.0). В связи с этим в сигнале между
участками с символами (ненулевыми) должны быть паузы с длительностью
минимум 7 мс, а сам полезный сигнал (для кодирования символа) должен
длиться минимум 80 мс.
__Файл с названием z100s.txt декодируется в *100#__
