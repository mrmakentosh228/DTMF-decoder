import matplotlib.pyplot as plt
import numpy as np
from math import cos, pi, sqrt


tabFreq = [697, 770, 852, 941, 1209, 1336, 1477, 1633]  # таблица частот
SymbolLenSec = 0.08  # длина буфера для сигнала одного символа в секундах
PauseLenSec = 0.007  # длина буфера для сигнала паузы в секундах
MinFreqAmplitude = 3.0  # минимальный порог для опознования амплитуды отдельной частоты
MaxPauseDev = 0.32  # максимальное стандартное отклонение сигнала в буфере паузы
                    # для опознования паузы
# словарь соответствий частот и кодируемых символов
symbols = {(697, 1209): '1',
           (697, 1336): '2',
           (697, 1477): '3',
           (697, 1633): 'A',
           (770, 1209): '4',
           (770, 1336): '5',
           (770, 1477): '6',
           (770, 1633): 'B',
           (852, 1209): '7',
           (852, 1336): '8',
           (852, 1477): '9',
           (852, 1633): 'C',
           (941, 1209): '*',
           (941, 1336): '0',
           (941, 1477): '#',
           (941, 1633): 'D'}


def goertzel(sample, sample_len, sample_freq, target_freq):  # реализация алгоритма Герцеля
    #  номер частоты, для которой вычисляется амплитуда в данной выборке
    k = int(target_freq / sample_freq * sample_len) + 1
    #  массив вычислеямых промежуточных значений для алгоритма
    s = np.zeros(sample_len + 1)
    #  просто постоянный коэффициент используемый в вычислениях
    super_k = 2 * pi * k / sample_len
    for i in range(2, sample_len + 1):
        s[i] = 2 * cos(super_k) * s[i - 1] - s[i - 2] + sample[i - 2]
    t1, t2 = s[-1], s[-2]
    amplitude = sqrt(t1 ** 2 - 2 * cos(super_k) * t1 * t2 + t2 ** 2)
    return amplitude


def input_processing(text):  # обработка входных данных
    f = open(text, mode='r')
    temp = f.readlines()
    num, freqx = map(int, temp[0].split())
    s = np.array(list(map(float, temp[1:])))
    f.close()
    #  возвращает массив отсчетов, их количество, и частоту дискретизации
    return s, num, freqx


def data_processing(data, num, freq):  # основная функция обработки массива данных
    #  длины буферов в отсчетах
    pause_len = int(PauseLenSec * freq)
    symbol_len = int(SymbolLenSec * freq)
    #  сами буфера
    pause = []
    symbol = []
    #  словарь, содержащий при обнаружении соответствующего символа счетчик его
    #  обнаружения
    found_symbols = {}
    #  конечный вывод распознанных символов
    output = []
    for i, s in enumerate(data):
        #  сначала заполняется буфер паузы
        if len(pause) < pause_len:
            pause.append(s)
        #  если находится в паузе (тишине) или на переходе с тишины на сигнал символа
        elif pause_check(pause):
            #  очищаются буфера, т. к. пауза найдена
            pause.clear()
            #  предыдущий блок отсчетов для детектирования символа уже был обработан либо не актуален
            symbol.clear()
            #  добавляем задетектированный символ в список для вывода
            make_symbol(found_symbols, output)
            #  не забыть добавить новый отсчет в буфер для паузы
            pause.append(s)
        else:  # не пауза
            #  перекидываем четверть старых отсчетов из буфера паузы в буфер символа
            while len(pause) > pause_len * 3 / 4 and len(symbol) < symbol_len:
                symbol.append(pause[0])
                pause.pop(0)
            #  если длины буфера символа достаточно для его обнаружения,
            if len(symbol) >= symbol_len:
                #  то вычисляем амплитуды для восьми необходимых частот
                amplitudes = find_amp(symbol, len(symbol), freq)
                #  находится наиболее возможный символ
                find_symb(amplitudes, found_symbols, i)
                #  удаляем четверть самых старых отсчетов
                while len(symbol) > symbol_len * 3 / 4:
                    symbol.pop(0)
            pause.append(s)
    return output


def make_symbol(comp_dict, comp_list):  # функция для добавления символа
    if len(comp_dict) > 0:
        #  символ с наибольшим счетчиком
        true_symbol = max(comp_dict, key=comp_dict.get)
        comp_list.append(true_symbol)
        #  очищаем словарь для поиска символов, т.к. сигнал, кодирующий символ завершен,
        #  либо только начинается
        comp_dict.clear()


def find_amp(symbol, num, freq):
    #  нахождение амплитуд восьми нужных частот в данной выборке
    #  (буфере символа)
    amplitudes = np.zeros(8)
    for i, f in enumerate(tabFreq):
        amplitudes[i] = goertzel(symbol, num, freq, f)
    return amplitudes


def find_symb(amplitudes, found_symbols, i):
    #  фунция находит две частоты с максимальными амплитудами и увеличивает
    #  счетчик символа, соответствующего этим амплитудам
    max1 = max(amplitudes)
    max1_ind = max(enumerate(amplitudes), key=lambda x: x[1])[0]
    amplitudes[max1_ind] = -1
    max2 = max(amplitudes)
    max2_ind = max(enumerate(amplitudes), key=lambda x: x[1])[0]
    amplitudes[max2_ind] = -1
    #  проверка на шум, если есть третья "большая" частота, то это сокрей всего он
    if max(amplitudes) > max2 / 2:
       return 1
    #  также проверяем амплитуды частот
    if max1 >= MinFreqAmplitude and max2 >= MinFreqAmplitude:
        freq1 = max(tabFreq[max1_ind], tabFreq[max2_ind])
        freq2 = min(tabFreq[max1_ind], tabFreq[max2_ind])
        symb = symbols[(freq2, freq1)]
        if symb not in found_symbols.keys():
            found_symbols[symb] = 0
        found_symbols[symb] += 1
        return freq2, freq1


def pause_check(pause):  # проверка буфера паузы для ее обнаружения
    npause = np.array(pause)
    #  считается стандартное отклонение отсчетов через дисперсию
    #  и сравнивается с пороговым значением
    deviation = sqrt(np.mean(npause*npause) - np.mean(npause) ** 2)
    return deviation <= MaxPauseDev

print("Введите название файла для расшифровки: ")
inp, N, frequence = input_processing(input())
out = data_processing(inp, N, frequence)
print("Результат декодировки: ")
print(''.join(out))


#  спектральные амплитуды всего сигнала для визуализации
amp = []
for l in tabFreq:
    amp.append(goertzel(inp, N, frequence, l))

plt.plot(tabFreq, amp, 'ro')
plt.vlines(tabFreq, np.zeros(8), amp)
plt.show()
