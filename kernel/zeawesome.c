/*
 *	Linux Magic System Request Key Hacks
 *
 *	(c) 1997 Martin Mares <mj@atrey.karlin.mff.cuni.cz>
 *	based on ideas by Pavel Machek <pavel@atrey.karlin.mff.cuni.cz>
 *
 *	(c) 2000 Crutcher Dunnavant <crutcher+kernel@datastacks.com>
 *	overhauled to use key registration
 *	based upon discusions in irc://irc.openprojects.net/#kernelnewbies
 *
 *	Copyright (c) 2010 Dmitry Torokhov
 *	Input handler conversion
 */

#define pr_fmt(fmt) KBUILD_MODNAME ": " fmt

#include <linux/sched.h>
#include <linux/interrupt.h>
#include <linux/mm.h>
#include <linux/fs.h>
#include <linux/mount.h>
#include <linux/kdev_t.h>
#include <linux/major.h>
#include <linux/reboot.h>

//#include <linux/sysrq.h>

#include <linux/kbd_kern.h>
#include <linux/proc_fs.h>
#include <linux/nmi.h>
#include <linux/quotaops.h>
#include <linux/perf_event.h>
#include <linux/kernel.h>
#include <linux/module.h>
#include <linux/suspend.h>
#include <linux/writeback.h>
#include <linux/buffer_head.h>		/* for fsync_bdev() */
#include <linux/swap.h>
#include <linux/spinlock.h>
#include <linux/vt_kern.h>
#include <linux/workqueue.h>
#include <linux/hrtimer.h>
#include <linux/oom.h>
#include <linux/slab.h>
#include <linux/input.h>
#include <linux/miscdevice.h>
#include <linux/semaphore.h>

#include <asm/ptrace.h>
#include <asm/irq_regs.h>

/* Whether we react on zeawesome keys or just ignore them */
static int __read_mostly zeawesome_enabled = true;
static bool __read_mostly zeawesome_always_enabled;

static bool zeawesome_on(void)
{
	return zeawesome_enabled || zeawesome_always_enabled;
}

/*
 * A value of 1 means 'all', other nonzero values are an op mask:
 */
static bool zeawesome_on_mask(int mask)
{
	return zeawesome_always_enabled ||
	       zeawesome_enabled == 1 ||
	       (zeawesome_enabled & mask);
}

static int __init zeawesome_always_enabled_setup(char *str)
{
	zeawesome_always_enabled = true;
	pr_info("zeawesome always enabled.\n");

	return 1;
}

__setup("zeawesome_always_enabled", zeawesome_always_enabled_setup);


static void zeawesome_handle_loglevel(int key)
{
	int i;

	i = key - '0';
	console_loglevel = 7;
	printk("Loglevel set to %d\n", i);
	console_loglevel = i;
}


/*
 * Signal zeawesome helper function.  Sends a signal to all user processes.
 */
static void send_sig_all(int sig)
{
	struct task_struct *p;

	for_each_process(p) {
		if (p->mm && !is_global_init(p))
			/* Not swapper, init nor kernel thread */
			force_sig(sig, p);
	}
}
/**
void handle_zeawesome(int key)
{
	if (zeawesome_on())
		__handle_zeawesome(key, true);
}
EXPORT_SYMBOL(handle_zeawesome);
*/
struct semaphore chrdev_sema;

static DEFINE_MUTEX(chrdev_mutex);
static char chrdev_buf[32];

static int zeawesome_chrdev_open(struct inode *inode, struct file *filp)
{
	if ((filp->f_mode & FMODE_READ) == 0)
		return -EINVAL;
	if (filp->f_mode & FMODE_WRITE)
		return -EINVAL;
	return 0;
}

static ssize_t zeawesome_chrdev_read(struct file *filp, char __user *buf,
			    size_t size, loff_t *offp)
{
//	ssize_t ret = 0;
	int err = 0;
//	int bytes_read, len;
	
	down_interruptible (&chrdev_sema);
	
	mutex_lock_interruptible(&chrdev_mutex);
	if (copy_to_user(buf, chrdev_buf, 1)) {
		err = -EFAULT;
		return 0;
	}
	memmove(chrdev_buf, chrdev_buf+1, sizeof(chrdev_buf)-1);
	mutex_unlock(&chrdev_mutex);
	return 1;
	
/*	while (size) {
		if (mutex_lock_interruptible(&rng_mutex)) {
			err = -ERESTARTSYS;
			goto out;
		}

		if (!current_rng) {
			err = -ENODEV;
			goto out_unlock;
		}

		if (!data_avail) {
			bytes_read = rng_get_data(current_rng, rng_buffer,
				sizeof(rng_buffer),
				!(filp->f_flags & O_NONBLOCK));
			if (bytes_read < 0) {
				err = bytes_read;
				goto out_unlock;
			}
			data_avail = bytes_read;
		}

		if (!data_avail) {
			if (filp->f_flags & O_NONBLOCK) {
				err = -EAGAIN;
				goto out_unlock;
			}
		} else {
			len = data_avail;
			if (len > size)
				len = size;

			data_avail -= len;

			if (copy_to_user(buf + ret, rng_buffer + data_avail,
								len)) {
				err = -EFAULT;
				goto out_unlock;
			}

			size -= len;
			ret += len;
		}

		mutex_unlock(&rng_mutex);

		if (need_resched())
			schedule_timeout_interruptible(1);

		if (signal_pending(current)) {
			err = -ERESTARTSYS;
			goto out;
		}
	}
out:
	return ret ? : err;
out_unlock:
	mutex_unlock(&rng_mutex);
	goto out;*/
	return 0;
}

static const struct file_operations zeawesome_chrdev_ops = {
	.owner	= THIS_MODULE,
	.open		= zeawesome_chrdev_open,
	.read		= zeawesome_chrdev_read,
	.llseek	= noop_llseek,
};

static struct miscdevice zeawesome_miscdev = {
	.minor	= MISC_DYNAMIC_MINOR,
	.name		= "zeawesome",
	.nodename	= "zeawesome",
	.fops		= &zeawesome_chrdev_ops,
};


#ifndef CONFIG_INPUT
#error CONFIG_INPUT required
#endif

static struct input_dev *mouse_dev;



/* Simple translation table for the SysRq keys */
static const unsigned char zeawesome_xlate[KEY_CNT] =
        "\000\0331234567890-=\177\t"                    /* 0x00 - 0x0f */
        "qwertyuiop[]\r\000as"                          /* 0x10 - 0x1f */
        "dfghjkl;'`\000\\zxcv"                          /* 0x20 - 0x2f */
        "bnm,./\000*\000 \000\201\202\203\204\205"      /* 0x30 - 0x3f */
        "\206\207\210\211\212\000\000789-456+1"         /* 0x40 - 0x4f */
        "230\177\000\000\213\214\000\000\000\000\000\000\000\000\000\000" /* 0x50 - 0x5f */
        "\r\000/";                                      /* 0x60 - 0x6f */


struct zeawesome_state {
	struct input_handle handle;
//	struct work_struct reinject_work;
	unsigned long key_down[BITS_TO_LONGS(KEY_CNT)];
	unsigned int alt;
	unsigned int alt_use;
	bool active;
	bool need_reinject;
	bool reinjecting;
};

/*
static void zeawesome_reinject_alt_zeawesome(struct work_struct *work)
{
	struct zeawesome_state *zeawesome =
			container_of(work, struct zeawesome_state, reinject_work);
	struct input_handle *handle = &zeawesome->handle;
	unsigned int alt_code = zeawesome->alt_use;

	if (zeawesome->need_reinject) {
		// we do not want the assignment to be reordered
		zeawesome->reinjecting = true;
		mb();

		// Simulate press and release of Alt + SysRq
		input_inject_event(handle, EV_KEY, alt_code, 1);
		input_inject_event(handle, EV_KEY, KEY_SYSRQ, 1);
		input_inject_event(handle, EV_SYN, SYN_REPORT, 1);

		input_inject_event(handle, EV_KEY, KEY_SYSRQ, 0);
		input_inject_event(handle, EV_KEY, alt_code, 0);
		input_inject_event(handle, EV_SYN, SYN_REPORT, 1);

		mb();
		zeawesome->reinjecting = false;
	}
}*/

static bool zeawesome_filter(struct input_handle *handle,
			 unsigned int type, unsigned int code, int value)
{
	struct zeawesome_state *zeawesome = handle->private;
	bool was_active = zeawesome->active;
	bool suppress;
	
	/*
	 * Do not filter anything if we are in the process of re-injecting
	 * Alt+SysRq combination.
	 */
	if (zeawesome->reinjecting)
		return false;

	switch (type) {

	case EV_SYN:
		suppress = false;
		break;

	case EV_KEY:
		if (code == KEY_CAPSLOCK) {
			suppress = true;
			if (value)
				zeawesome->active = true;
			else
				zeawesome->active = false;
			break;
		}else if (zeawesome->active) {
			suppress = true;
			switch (code) {
			case KEY_F:		code = BTN_LEFT;		goto rep_key;
			case KEY_D:		code = BTN_RIGHT;		goto rep_key;
			case KEY_S:		code = BTN_MIDDLE;	goto rep_key;
			
			case KEY_C:		code = BTN_SIDE;		goto rep_key;
			case KEY_V:		code = BTN_EXTRA;		goto rep_key;
			rep_key:
				if (value == 2)
					break;
				input_report_key(mouse_dev, code, value);
				input_sync(mouse_dev);
				break;
			
			case KEY_X:
				if (value == 0)
					break;
				input_report_rel(mouse_dev, REL_WHEEL, -1);
				input_sync(mouse_dev);
				break;
			case KEY_E:
				if (value == 0)
					break;
				input_report_rel(mouse_dev, REL_WHEEL, 1);
				input_sync(mouse_dev);
				break;
			default:
				if (value == 2)
					break;
				if (code >= KEY_CNT)
					break;
				if (	code == KEY_LEFTALT || code == KEY_RIGHTALT
					|| code == KEY_LEFTCTRL || code == KEY_RIGHTCTRL
					|| code == KEY_LEFTSHIFT || code == KEY_RIGHTSHIFT
					|| code == KEY_LEFTMETA || code == KEY_RIGHTMETA
				/*	|| ( !(	(code >= KEY_1 && code <= KEY_EQUAL)
							&& (code >= KEY_TAB && code <= KEY_P)
							&& (code >= KEY_TAB && code <= KEY_P)
					))/**/
				) {
					suppress = false;
					break;
				}
				mutex_lock_interruptible(&chrdev_mutex);
				if (chrdev_sema.count < 32) {
				//	chrdev_buf[chrdev_sema.count] = zeawesome_xlate[code] | (value << 7);
					chrdev_buf[chrdev_sema.count] = code | (value << 7);
					mutex_unlock(&chrdev_mutex);
					up(&chrdev_sema);
				}else {
					mutex_unlock(&chrdev_mutex);
				}
				break;
			}
			break;
		}else {
			suppress = false;
			break;
		}
		
		switch (code) {

		case KEY_LEFTALT:
		case KEY_RIGHTALT:
			if (!value) {
				/* One of ALTs is being released */
				if (zeawesome->active && code == zeawesome->alt_use)
					zeawesome->active = false;

				zeawesome->alt = KEY_RESERVED;

			} else if (value != 2) {
				zeawesome->alt = code;
				zeawesome->need_reinject = false;
			}
			break;

		case KEY_SYSRQ:
			if (value == 1 && zeawesome->alt != KEY_RESERVED) {
				zeawesome->active = true;
				zeawesome->alt_use = zeawesome->alt;
				/*
				 * If nothing else will be pressed we'll need
				 * to re-inject Alt-SysRq keysroke.
				 */
				zeawesome->need_reinject = true;
			}

			/*
			 * Pretend that zeawesome was never pressed at all. This
			 * is needed to properly handle KGDB which will try
			 * to release all keys after exiting debugger. If we
			 * do not clear key bit it KGDB will end up sending
			 * release events for Alt and SysRq, potentially
			 * triggering print screen function.
			 */
			if (zeawesome->active)
				clear_bit(KEY_SYSRQ, handle->dev->key);

			break;

		default:
			if (zeawesome->active && value && value != 2) {
				zeawesome->need_reinject = false;
			//	__handle_zeawesome(zeawesome_xlate[code], true);
			}
			break;
		}

		suppress = zeawesome->active;

		if (!zeawesome->active) {
			/*
			 * If we are not suppressing key presses keep track of
			 * keyboard state so we can release keys that have been
			 * pressed before entering SysRq mode.
			 */
			if (value)
				set_bit(code, zeawesome->key_down);
			else
				clear_bit(code, zeawesome->key_down);
			
		//	if (was_active)
		//		schedule_work(&zeawesome->reinject_work);
			
		} else if (value == 0 &&
			   test_and_clear_bit(code, zeawesome->key_down)) {
			/*
			 * Pass on release events for keys that was pressed before
			 * entering SysRq mode.
			 */
			suppress = false;
		}
		break;

	default:
		suppress = zeawesome->active;
		break;
	}

	return suppress;
}

static int zeawesome_connect(struct input_handler *handler,
			 struct input_dev *dev,
			 const struct input_device_id *id)
{
	struct zeawesome_state *zeawesome;
	int error;
	printk(KERN_INFO "zeawesome_connect\n");
	
	zeawesome = kzalloc(sizeof(struct zeawesome_state), GFP_KERNEL);
	if (!zeawesome)
		return -ENOMEM;
	
//	INIT_WORK(&zeawesome->reinject_work, zeawesome_reinject_alt_zeawesome);
	
	zeawesome->handle.dev = dev;
	zeawesome->handle.handler = handler;
	zeawesome->handle.name = "zeawesome";
	zeawesome->handle.private = zeawesome;

	error = input_register_handle(&zeawesome->handle);
	if (error) {
		pr_err("Failed to register input zeawesome handler, error %d\n",
			error);
		goto err_free;
	}

	error = input_open_device(&zeawesome->handle);
	if (error) {
		pr_err("Failed to open input device, error %d\n", error);
		goto err_unregister;
	}

	return 0;

 err_unregister:
	input_unregister_handle(&zeawesome->handle);
 err_free:
	kfree(zeawesome);
	return error;
}

static void zeawesome_disconnect(struct input_handle *handle)
{
	struct zeawesome_state *zeawesome = handle->private;
	
	input_close_device(handle);
//	cancel_work_sync(&zeawesome->reinject_work);
	input_unregister_handle(handle);
	kfree(zeawesome);
}

/*
 * We are matching on KEY_LEFTALT instead of KEY_SYSRQ because not all
 * keyboards have SysRq key predefined and so user may add it to keymap
 * later, but we expect all such keyboards to have left alt.
 */
static const struct input_device_id zeawesome_ids[] = {
	{
		.flags = INPUT_DEVICE_ID_MATCH_EVBIT |
				INPUT_DEVICE_ID_MATCH_KEYBIT,
		.evbit = { BIT_MASK(EV_KEY) },
		.keybit = { BIT_MASK(KEY_CAPSLOCK) },
	},
	{ },
};

static struct input_handler zeawesome_handler = {
	.filter	= zeawesome_filter,
	.connect	= zeawesome_connect,
	.disconnect	= zeawesome_disconnect,
	.name		= "zeawesome",
	.id_table	= zeawesome_ids,
};

static bool zeawesome_handler_registered;

static inline void zeawesome_register_handler(void)
{
	int error;
	printk(KERN_INFO "zeawesome_register_handler\n");
	
	error = input_register_handler(&zeawesome_handler);
	if (error)
		pr_err("Failed to register input handler, error %d", error);
	else
		zeawesome_handler_registered = true;
}

static inline void zeawesome_unregister_handler(void)
{
	printk(KERN_INFO "zeawesome_unregister_handler\n");
	if (zeawesome_handler_registered) {
		input_unregister_handler(&zeawesome_handler);
		zeawesome_handler_registered = false;
	}
}

//#else
/*
static inline void zeawesome_register_handler(void)
{
}

static inline void zeawesome_unregister_handler(void)
{
}
*/

//#endif /* CONFIG_INPUT */

int zeawesome_toggle_support(int enable_mask)
{
	bool was_enabled = zeawesome_on();

	zeawesome_enabled = enable_mask;

	if (was_enabled != zeawesome_on()) {
		if (zeawesome_on())
			zeawesome_register_handler();
		else
			zeawesome_unregister_handler();
	}

	return 0;
}

#ifdef CONFIG_PROC_FS
/*
 * writing 'C' to /proc/zeawesome-trigger is like zeawesome-C
 */
static ssize_t write_zeawesome_trigger(struct file *file, const char __user *buf,
				   size_t count, loff_t *ppos)
{
	if (count) {
		char c;

		if (get_user(c, buf))
			return -EFAULT;
		if (c == 'k') {
			zeawesome_unregister_handler();
		}
	//	__handle_zeawesome(c, false);
	}

	return count;
}

static const struct file_operations proc_zeawesome_trigger_operations = {
	.write		= write_zeawesome_trigger,
	.llseek		= noop_llseek,
};

static void zeawesome_init_procfs(void)
{
	if (!proc_create("zeawesome-trigger", S_IWUSR, NULL,
			 &proc_zeawesome_trigger_operations))
		pr_err("Failed to register proc interface\n");
}

static void zeawesome_exit_procfs(void)
{
	remove_proc_entry("zeawesome-trigger", NULL);
}

#else

static inline void zeawesome_init_procfs(void)
{
}

#endif /* CONFIG_PROC_FS */

static int __init zeawesome_init(void)
{
	int error;
	
	printk(KERN_INFO "zeawesome_init\n");
	
	zeawesome_init_procfs();
	
	sema_init (&chrdev_sema, 0);
	
//	if (request_irq(BUTTON_IRQ, button_interrupt, 0, "button", NULL)) {
//		printk(KERN_ERR "button.c: Can't allocate irq %d\n", button_irq);
//		return -EBUSY;
//	}
	
	error = misc_register(&zeawesome_miscdev);
	if (error)
		return error;
	
	printk("zeawesome: got minor %i\n", zeawesome_miscdev.minor);
	
	
	mouse_dev = input_allocate_device();
	if (!mouse_dev) {
		printk(KERN_ERR "zeawesome: Not enough memory\n");
		error = -ENOMEM;
		goto err_free_irq;
	}

	mouse_dev->evbit[0] = BIT_MASK(EV_KEY) |  BIT_MASK(EV_REL);
	
	mouse_dev->relbit[BIT_WORD(REL_X)] |= BIT_MASK(REL_X);
	mouse_dev->relbit[BIT_WORD(REL_Y)] |= BIT_MASK(REL_Y);
	mouse_dev->relbit[BIT_WORD(REL_WHEEL)] |= BIT_MASK(REL_WHEEL);
	
//	mouse_dev->keybit[BIT_WORD(BTN_0)] = BIT_MASK(BTN_0);
	#define dkey(key) mouse_dev->keybit[BIT_WORD(key)] |= BIT_MASK(key)
	dkey(BTN_LEFT);
	dkey(BTN_RIGHT);
	dkey(BTN_MIDDLE);
	
	dkey(BTN_SIDE);
	dkey(BTN_EXTRA);
	
	#undef dkey
	mouse_dev->name = "Generic Mouse";
	
	error = input_register_device(mouse_dev);
	if (error) {
		printk(KERN_ERR "zeawesome: Failed to register device\n");
		goto err_free_dev;
	}
	
	if (zeawesome_on())
		zeawesome_register_handler();
	
	return 0;
	
err_free_dev:
	input_free_device(mouse_dev);
err_free_irq:
//	free_irq(BUTTON_IRQ, button_interrupt);
	return error;
}
module_init(zeawesome_init);

static void __exit zeawesome_exit(void)
{
	zeawesome_exit_procfs ();
	
	misc_deregister(&zeawesome_miscdev);
	
	input_unregister_device(mouse_dev);
	
	zeawesome_unregister_handler();
}

module_exit(zeawesome_exit);
