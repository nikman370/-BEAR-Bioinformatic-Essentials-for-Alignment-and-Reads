def check_gc(sequence: str, gc_bounds: tuple[float, float] | float) -> bool:
    """
    Проверяет: укладывается ли GC-состав последовательности в заданные рамки

    Arguments:
        sequence: Строка нуклеотидов (ATGC)
        gc_bounds: Кортеж (мин, макс) ИЛИ одно число (макс)

    Returns:
        True, если рид проходит фильтр, иначе False.
    """

    seq_upper = sequence.upper()  # Переводим G и C (в верхнем регистр)
    gc_count = seq_upper.count("G") + seq_upper.count("C")  # Считаем сколько GC

    # вдруг нам не введут рид, чтобы не делить на ноль, если рид пустой
    if not sequence:
        gc_pct = 0.0
    else:
        gc_pct = (gc_count / len(sequence)) * 100

    # проверяем, какой тип данных нам передали в gc_bounds
    if type(gc_bounds) is tuple:
        min_gc = gc_bounds[0]  # нижняя граница — первое число в кортеже
        max_gc = gc_bounds[1]  # верхняя граница — второе число в кортеже
    else:
        min_gc = 0.0  # нижней границы нет, значит это 0%
        max_gc = gc_bounds  # само число становится верхней границей

    return (
        min_gc <= gc_pct <= max_gc
    )  # истина - если число внутри границ. или False - если число вылетело за границы


def check_length(sequence: str, length_bounds: tuple[int, int] | int) -> bool:
    # Проверяем, укладывается ли длина последовательности в заданные границы.
    seq_len = len(sequence)

    if type(length_bounds) is tuple:
        min_len = length_bounds[0]
        max_len = length_bounds[1]
    else:
        min_len = 0
        max_len = length_bounds

    return min_len <= seq_len <= max_len


def check_quality(quality_str: str, quality_threshold: float) -> bool:
    """
    Проверяет, превышает ли среднее качество рида заданный порог Phred33.

    Arguments:
        quality_str: Строка качества из FASTQ файла.
        quality_threshold: Минимально допустимое среднее качество.

    Returns:
        True, если среднее качество >= порога, иначе False.
    """
    if not quality_str:
        return False

    total_quality = 0  # переменная для суммы

    for char in quality_str:  # Перебираем каждый символ в строке качества
        phred_score = ord(char) - 33  # Переводим символ в число Phred33
        total_quality = total_quality + phred_score  # Добавляем в сумму

    mean_quality = total_quality / len(quality_str)  # Считаем среднее

    return mean_quality >= quality_threshold  # Проверяем, проходит ли по порогу
